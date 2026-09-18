#!/usr/bin/env python3
# tiny 2-leg spinal CPG, reduced network: k=0.25 connectivity, m=0.40 neurons.
# STDP-CAPABLE: adds optional plasticity on cut/bs->RG inputs (--stdp-lambda).
# Reciprocal inhibition stays static, so the flexor-sustaining alternation is preserved.
# --stdp-lambda 0 reproduces cpg_final_k025_m040.py exactly.
# Defaults below ARE the validated recipe -- no flags needed to reproduce it.
#   inh_comp=3.0  X0_E=1.10  K_E=15  X0_F=1.42  K_F=30
# Result: corr(Force-E,Force-F) = -0.911 +/- 0.073 (n=8 seeds), baseline k=1.00 = -0.953.
# Synapses ~= m^2*k = 4% of baseline; neurons = 40%.
# Run with --threads 1: NEST is nondeterministic at higher thread counts for these sizes.
# -*- coding: utf-8 -*-
"""
cpg_variantQ_k050_m040_fixedweights.py  —  VARIANT Q
=====================================================
k=0.50, m=0.40, STDP OFF.

THE I -> P STORY, IN ONE LINE:
  Every variant since I improved corr(RG) (phase alignment) without
  ever moving Force-F's trough. N validated the RG circuit itself
  against the real baseline (I_E_RGF=0.9, W_INE2RGF=-8.0,
  W_RG_REC_F=5.5, paced-gait filters 40/80) and got corr(Force) into
  the neighborhood of baseline (-0.8835 mean) -- but the Force-F/
  Activation-F plots still showed a floor stuck near ceiling. O
  (motor-pool asymmetry) and P (RG_RECIP_E density, re-tested on N's
  corrected foundation) both failed to move it (-0.8629, -0.8575).
  Three different layers, three failures, same visual symptom. That
  consistency is the actual finding: the bug isn't in the RG circuit,
  the motor pool, or the connectivity density. It's downstream of all
  of them, in the activation gate itself.

THE BUG:
    de = 1/(1+exp(-8*(r_rge/100 - 0.6)))
    df = 1/(1+exp(-8*(r_rgf/100 - 0.6)))
  ACT_GATE_X0=0.6 is SHARED between E and F, inherited symmetrically
  from a model where nothing else is symmetric: F->E inhibition is
  2x-6x stronger than E->F (by design, per baseline's own Zhang-ratio
  docstring), F's drive is deliberately lower than E's (I_E_RGF=0.9 <
  I_E_RGE=1.0, also by design), and F's neuron population sees fewer
  INF interneurons at m=0.40. There is no reason RG-F's ABSOLUTE
  firing-rate range should coincide with RG-E's just because their
  PHASES are locked in antiphase (corr(RG)~-0.97 across I-P proves
  the phase lock is real and has been real since I). Correlation is
  scale/threshold-invariant; a hard-coded gate crossing point is not.
  A population that phase-locks correctly but simply operates at a
  lower absolute rate range will phase-lock perfectly on corr(RG) and
  STILL never cross a threshold calibrated to the other population's
  range -- which is exactly the "RG-F oscillates, force never drops"
  pattern seen in every plot since I.

Variant Q strategy: stop assuming E and F share a gate threshold just
because they share a formula. Split ACT_GATE_X0 into independent
per-population values:

  ACT_GATE_X0_E : 0.6   (unchanged -- E's range is already correctly
                          gated; corr(Force) contributions from the E
                          side have looked right throughout I-P)
  ACT_GATE_X0_F : 0.3   (NEW -- lower threshold, calibrated to F's
                          actual lower operating range instead of
                          inheriting E's. This is a genuinely new
                          lever: nothing since Variant H has touched
                          the gate itself, only what feeds into r_rgf.)

Everything else held at N's validated values (the last confirmed-good
RG-circuit state; O's motor-pool split and P's density re-test are
both dropped since neither helped):
  I_E_RGF=0.9, W_INE2RGF=-8.0, W_INF2RGE=-48.0, W_RG_REC_F=5.5,
  W_MOTOR_RECIP=-22.0 (single symmetric value, O's split dropped),
  P_RG_RECIP_E=0.15*_k, P_RG_RECIP_F=0.30*_k (P's density bump
  dropped), TAU_ACT/TAU_FORCE=40/40, 80/80 (paced-gait defaults,
  confirmed correct against baseline in N), Ia-F ablated
  (IA_BASE_HZ_F=5.0, IA_K_FORCE_F=IA_K_STRETCH_F=0.0).

Expected result:
  corr(RG) should stay ~-0.96 to -0.97 -- this change touches only the
  gate's crossing point, not the spiking circuit that produces phase.
  Activation-F and Force-F should show a REAL trough now, closing
  toward E's range whenever r_rgf dips below 30 Hz (X0=0.3 * ref=100)
  instead of needing 60 Hz (X0=0.6) that it may rarely reach.
  corr(Force) target: should move for the first time since N,
  ideally climbing toward/past baseline's -0.90.
  WATCH FOR: if ACT_GATE_X0_F=0.3 is too permissive, the gate could
  stay OPEN too much of the cycle and Force-F could fail to rise
  properly during F's active phase (a new failure mode, opposite of
  the one we've been chasing) -- check that Force-F still reaches a
  comparable PEAK to Force-E, not just a lower trough. If the peak
  degrades, the next step is a value between 0.3 and 0.6 (e.g. 0.45)
  rather than pushing lower.
"""
import argparse, os, time
from datetime import datetime
import numpy as np, h5py, nest

LEGS = ("L", "R")

# ── Scale factors ──────────────────────────────────────────
_m, _k = 0.40, 0.25   # VARIANT T: k=0.25, m=0.40 (deeper connectivity cut, same gate+IaF fixes as S)

N_CUT    = int(100*_m)    # 40
N_BS     = int(100*_m)    # 40
N_RG_TOTAL = int(200*_m)  # 80
N_RG_E   = N_RG_TOTAL//2  # 40
N_RG_F   = N_RG_TOTAL - N_RG_E  # 40
N_MOTOR_E= int(100*_m)    # 40
N_MOTOR_F= int(100*_m)    # 40
N_MUS_E  = int(100*_m)    # 40
N_MUS_F  = int(100*_m)    # 40
N_IA_E   = int(100*_m)    # 40
N_IA_F   = int(100*_m)    # 40
N_IA_INT = int(50*_m)     # 20
N_INE    = int(50*_m)     # 20
N_INF    = int(50*_m)     # 20

# ── Weights (N's validated RG circuit; O/P deltas dropped) ─
W_IA_IN2INT   =  6.0
W_IA_INT2ANT  = -10.0
W_IA2IN       =  6.0
W_FLEX_AFF2RGF=  6.0
W_RG2INE      =  12.0;  W_RG2INF  =  18.0
W_INE2RGF     =  -8.0   # N's reverted baseline-matched value (Zhang 6:1 ratio)
W_INF2RGE     = -48.0
W_RG_REC_E    =   4.0;  W_RG_REC_F=   5.5   # N: reverted to baseline's 5.5
W_MOTOR_RECIP = -22.0   # N/H value; O's directional split dropped (no benefit)
W_M2MUS       =   1.0
W_CUT2INE     =   6.0
W_COMM_F_INH  = -20.0;  W_COMM_E_INH = -8.0
W0_IN = 22.0
W0_RM = 30.0
BASE_DRIVE_W  = 1.0

# ── Connectivity (× k=0.50) ────────────────────────────────
P_IN_STDP      = 0.5
P_RG_REC       = 0.12*_k      # 0.06
P_RG_RECIP_F   = 0.30*_k      # 0.15
P_RG_RECIP_E   = 0.15*_k      # 0.075 -- N's value; P's 0.30*_k density bump dropped
P_MOTOR_RECIP  = 0.25*_k      # 0.125
P_M2MUS        = 0.80*_k      # 0.40
P_IA2IN        = 0.25*_k      # 0.125
P_COMM_F       = 0.22*_k      # 0.11
P_COMM_E       = 0.10*_k      # 0.05
P_CUT2INE      = 0.30
P_CUT2RGE_STATIC = 0.35
IA2RG_P        = 0.4
BASE_DRIVE_P   = 0.10

# ── BS / Ia drive ──────────────────────────────────────────
N_PHASES          = 6
CUT_RATE_ON_HZ    = 100.0; CUT_RATE_OFF_HZ = 0.0
BS_REGULAR_HZ     = 60.0
BS_REGULAR_DESYNC = "linspace"; BS_REGULAR_JITTER_MS = 0.3
BS_RATE_MIN_HZ    = 0.0; BS_RATE_BASE_HZ = 0.0; BS_NOISE_STD_HZ = 0.0
LEFT_RIGHT_BIAS_IE = 0.12
ENABLE_COMMISSURAL = True

# ── Delays ─────────────────────────────────────────────────
DELAY_MS = 1.0; DELAY_RECIP_MS = 1.0
DELAY_MOTOR_RECIP_E2F_MS = 1.5; DELAY_MOTOR_RECIP_F2E_MS = 1.0
DELAY_COMM_MS = 1.0
DELAY_PRESETS = {
    "rat": {
        "cut_to_rg":{"syn_delay_ms":1.0,"length_m":0.005,"velocity_mps":5.0},
        "bs_to_rg": {"syn_delay_ms":1.0,"length_m":0.010,"velocity_mps":5.0},
        "base_to_rg":{"syn_delay_ms":1.0,"length_m":0.010,"velocity_mps":5.0},
        "rg_to_m":  {"syn_delay_ms":1.0,"length_m":0.005,"velocity_mps":5.0},
        "m_to_mus": {"syn_delay_ms":1.0,"length_m":0.005,"velocity_mps":5.0},
        "ia_path":  {"syn_delay_ms":1.0,"length_m":0.010,"velocity_mps":10.0},
        "rg_rec":   {"syn_delay_ms":0.8,"length_m":0.002,"velocity_mps":5.0},
        "rg_recip": {"syn_delay_ms":1.0,"length_m":0.004,"velocity_mps":5.0},
        "motor_e2f":{"syn_delay_ms":1.0,"length_m":0.004,"velocity_mps":5.0},
        "motor_f2e":{"syn_delay_ms":1.0,"length_m":0.004,"velocity_mps":5.0},
        "commissural":{"syn_delay_ms":1.0,"length_m":0.006,"velocity_mps":5.0},
    }
}

def _mk_delay(model, species, key, fallback, res, jitter, scale):
    if model == "fixed": return float(fallback)
    p = DELAY_PRESETS.get(species, DELAY_PRESETS["rat"])
    pr= p.get(key, None)
    if pr is None: base = float(fallback)
    else: base = scale*(float(pr["syn_delay_ms"])+(float(pr["length_m"])/max(1e-9,float(pr["velocity_mps"]))*1000))
    d = float(base)
    if jitter > 0: d = d + nest.random.normal(mean=0.0, std=jitter)
    try:    d = nest.math.max(d, float(max(1e-9, res)))
    except: d = float(max(float(d), float(max(1e-9, res))))
    return d

# ── Izhikevich / muscle params ──────────────────────────────
izh_params     = dict(a=0.02,b=0.2,c=-65.0,d=8.0,V_th=30.0,V_min=-120.0)
izh_inh_params = dict(a=0.1, b=0.2,c=-65.0,d=2.0,V_th=30.0,V_min=-120.0)
I_E_RGE=1.0; I_E_RGF=0.9; I_E_MOTOR=1.0   # N: RGF drive reverted to baseline's 0.9
RGF_A=0.02; RGF_B=0.2; RGF_C=-55.0; RGF_D=4.0
TAU_ACT_RISE_MS=20.0; TAU_ACT_DECAY_MS=20.0
ACT_MAX=1.2; ACT_SAT_K=0.02

# ── VARIANT Q: THE NEW LEVER ────────────────────────────────
# ACT_GATE_X0 was a single shared threshold (0.6) applied identically to
# both r_rge/100 and r_rgf/100 in the activation-gate sigmoid. Nothing
# else about E vs F is symmetric (drive, inhibition ratio, population),
# so there is no reason their absolute firing-rate ranges should share
# a crossing point just because their phases are locked (corr(RG)~-0.97
# has been true since Variant I without ever moving Force-F). Splitting
# this is a genuinely new lever -- no prior variant touched the gate.
ACT_GATE_K     = 8.0   # steepness for E's gate (unchanged -- E's ~40-230Hz range is wide)
ACT_GATE_K_F   = 15.0  # NEW -- steeper slope for F's gate, since RG-F's dynamic range
                        # (~110-150Hz, measured from rg_rate plots) is much narrower than
                        # E's (~40-230Hz); K=8 alone would leave F stuck in the sigmoid's
                        # shallow middle region even with a correctly centered threshold.
ACT_GATE_X0_E  = 0.6   # unchanged -- E's range already gates correctly
ACT_GATE_X0_F  = 1.30  # CORRECTED (was 0.3) -- 0.3 assumed RG-F could dip near 30Hz, but
                        # rg_rate diagnostic plots show RG-F actually operates at ~110-150Hz
                        # in this k=0.50/m=0.40 network (never bursts low like RG-E does),
                        # so the old threshold sat ~80-120Hz BELOW the entire trace, leaving
                        # the gate permanently saturated open (Activation-F pinned ~1.1-1.2).
                        # 1.30 = midpoint of the measured 110-150Hz range, normalized /100.

TAU_FORCE_RISE_MS=30.0; TAU_FORCE_DECAY_MS=30.0
FORCE_MAX=25.0; FORCE_SAT_K=1.0
TAU_LENGTH_MS=260.0; L0=1.0; L_MIN=0.5; L_MAX=2.0
SHORTEN_GAIN=0.010; STRETCH_GAIN=0.35
IA_BASE_HZ_E=10.0; IA_K_FORCE_E=6.0; IA_K_STRETCH_E=250.0
IA_BASE_HZ_F=5.0;  IA_K_FORCE_F=6.0; IA_K_STRETCH_F=250.0
# RESTORED (was IA_K_FORCE_F=IA_K_STRETCH_F=0.0, the "IaF_ablated" state carried over
# from Variant N/Q/R): Ia-F was hardcoded to a flat 5Hz regardless of flexor force or
# stretch, completely disconnecting the flexor's proprioceptive feedback loop. Mirrored
# to E's values here as the starting point -- E's gains were already validated correct
# against baseline. Override with --ia-k-force-f / --ia-k-stretch-f to retune F
# independently if 6.0/250.0 over- or under-shoots F's (different) force/length range.
IA_RATE_MAX_HZ=500.0

# gate outputs, logged so we can check each gate actually opens and closes
_GATE_DIAG = {"de": [], "df": []}

def clamp(x,lo,hi): return float(max(lo,min(hi,x)))
def new_spikes(rec,last):
    cur=int(nest.GetStatus(rec,"n_events")[0]); return cur-last, cur
def safe_len(model):
    try: return len(nest.GetConnections(synapse_model=model))
    except: return -1

# ═══════════════════════════════════════════════════════════
def main():
    global TAU_ACT_RISE_MS,TAU_ACT_DECAY_MS,TAU_FORCE_RISE_MS,TAU_FORCE_DECAY_MS
    global W_IA2IN,W_INF2RGE,W_INE2RGF,W_COMM_F_INH,W_COMM_E_INH,BS_RATE_BASE_HZ
    ap=argparse.ArgumentParser(description="CPG fixed-weights variant Q, k=0.50 m=0.40")
    ap.add_argument("--out",              default="cpg_run.h5")
    ap.add_argument("--outdir",           default=".")
    ap.add_argument("--tag",              default="cpg_varQ")
    ap.add_argument("--run-name",         default="")
    ap.add_argument("--seed",    type=int,default=12345)
    ap.add_argument("--sweep-pairs",      default="")
    ap.add_argument("--sweep-run-idx",type=int,default=-1)
    ap.add_argument("--sweep-dist",       default="lognormal_cv")
    ap.add_argument("--sim-ms",  type=float,default=30000.0)
    ap.add_argument("--dt-ms",   type=float,default=10.0)
    ap.add_argument("--threads", type=int,  default=10)
    ap.add_argument("--print-every",type=int,default=50)
    ap.add_argument("--weight-sample-ms",type=float,default=100.0)
    ap.add_argument("--rate-update-ms",  type=float,default=20.0)
    ap.add_argument("--resolution-ms",   type=float,default=0.2)
    ap.add_argument("--simulate-chunk-ms",type=float,default=50.0)
    ap.add_argument("--long-run",action="store_true")
    ap.add_argument("--species",          default="rat",choices=["rat","human"])
    ap.add_argument("--delay-model",      default="fixed",choices=["fixed","length_velocity"])
    ap.add_argument("--delay-jitter-ms",type=float,default=0.2)
    ap.add_argument("--delay-scale",    type=float,default=1.0)
    ap.add_argument("--max-weight-conns",type=int, default=0)
    ap.add_argument("--save-weights",     default="none",choices=["none","final","snapshots"])
    ap.add_argument("--nest-verbosity",   default="M_ERROR")
    ap.add_argument("--bs-base-hz",     type=float,default=0.0)
    ap.add_argument("--bs-noise-std-hz",type=float,default=0.0)
    ap.add_argument("--enforce-tonic-bs",action="store_true")
    ap.add_argument("--paced-gait",      action="store_true")
    ap.add_argument("--step-period-ms",type=float,default=520.0)
    ap.add_argument("--stance-fraction",type=float,default=0.5)
    ap.add_argument("--n-ia-groups",   type=int,  default=3)
    ap.add_argument("--ia-ext-hz",     type=float,nargs="+",default=[60.,80.,100.])
    ap.add_argument("--ia-ext-f-hz",   type=float,default=80.0)
    ap.add_argument("--ia-feedback-gain",type=float,default=1.0)
    ap.add_argument("--cut-feedback-gain",type=float,default=1.0)
    ap.add_argument("--static-weight-cv",type=float,default=0.5)
    ap.add_argument("--ablate-ia-loop", action="store_true")
    ap.add_argument("--ablate-asym",    action="store_true")
    ap.add_argument("--ablate-comm",    action="store_true")
    ap.add_argument("--act-gate-x0-f", type=float, default=1.42)  # F gate threshold
    ap.add_argument("--act-gate-k-f", type=float, default=30.0)   # F gate steepness
    ap.add_argument("--ia-k-force-f",   type=float, default=None,
                     help="Override IA_K_FORCE_F (restored Ia-F force-feedback gain, default 6.0)")
    ap.add_argument("--stdp-lambda", type=float, default=0.0,
                    help="STDP learning rate on the sensory/descending inputs to RG. "
                         "0 = static weights (identical to the validated model). >0 makes "
                         "cut->RG-E and bs->RG plastic via stdp_synapse. The reciprocal "
                         "inhibition that sustains alternation is NEVER made plastic.")
    ap.add_argument("--stdp-wmax", type=float, default=40.0,
                    help="Upper weight bound for plastic input synapses (pA)")
    ap.add_argument("--k-conn", type=float, default=0.25)      # synaptic connectivity scale
    ap.add_argument("--m-neurons", type=float, default=0.40)   # neuron-count scale
    ap.add_argument("--act-gate-x0-e", type=float, default=1.10)  # E gate threshold
    ap.add_argument("--act-gate-k-e", type=float, default=15.0)   # E gate steepness
    ap.add_argument("--inh-comp", type=float, default=3.0)     # E->F inhibition gain, holds p*w constant
    ap.add_argument("--inh-comp-f", type=float, default=1.0)   # F->E inhibition gain (uncompensated)
    ap.add_argument("--w-ine2rgf", type=float, default=None)   # direct weight override
    ap.add_argument("--w-inf2rge", type=float, default=None)   # direct weight override
    ap.add_argument("--i-e-rgf", type=float, default=None)     # RG-F tonic drive
    ap.add_argument("--w-rg-rec-f", type=float, default=None)  # RG-F recurrent excitation
    ap.add_argument("--ia-k-stretch-f", type=float, default=None,
                     help="Override IA_K_STRETCH_F (restored Ia-F stretch-feedback gain, default 250.0)")
    args=ap.parse_args()

    # k/m are module-level and frozen at import: every derived constant must be recomputed here
    global _k, _m
    global N_CUT, N_BS, N_RG_TOTAL, N_RG_E, N_RG_F, N_MOTOR_E, N_MOTOR_F
    global N_MUS_E, N_MUS_F, N_IA_E, N_IA_F, N_IA_INT, N_INE, N_INF
    global P_RG_REC, P_RG_RECIP_F, P_RG_RECIP_E, P_MOTOR_RECIP
    global P_M2MUS, P_IA2IN, P_COMM_F, P_COMM_E
    _k = float(args.k_conn); _m = float(args.m_neurons)
    if True:
        N_CUT = N_BS = int(100*_m)
        N_RG_TOTAL = int(200*_m)
        N_RG_E = N_RG_TOTAL//2
        N_RG_F = N_RG_TOTAL - N_RG_E
        N_MOTOR_E = N_MOTOR_F = int(100*_m)
        N_MUS_E   = N_MUS_F   = int(100*_m)
        N_IA_E    = N_IA_F    = int(100*_m)
        N_IA_INT  = N_INE = N_INF = int(50*_m)
    if True:
        P_RG_REC      = 0.12*_k
        P_RG_RECIP_F  = 0.30*_k
        P_RG_RECIP_E  = 0.15*_k
        P_MOTOR_RECIP = 0.25*_k
        P_M2MUS       = 0.80*_k
        P_IA2IN       = 0.25*_k
        P_COMM_F      = 0.22*_k
        P_COMM_E      = 0.10*_k
    print(f"[varAB] k={_k} m={_m}  N_RG_E={N_RG_E} N_INE={N_INE} N_MOTOR_E={N_MOTOR_E}"
          f"  P_RG_RECIP_E={P_RG_RECIP_E:.4f} P_M2MUS={P_M2MUS:.4f}")
    global ACT_GATE_X0_F, ACT_GATE_K_F, IA_K_FORCE_F, IA_K_STRETCH_F
    global ACT_GATE_X0_E, ACT_GATE_K
    global W_RG2INE, W_RG2INF, W_RG_REC_F, I_E_RGF
    if args.act_gate_x0_e is not None:
        ACT_GATE_X0_E = float(args.act_gate_x0_e)
    if args.act_gate_k_e is not None:
        ACT_GATE_K = float(args.act_gate_k_e)
    if args.act_gate_x0_f is not None:
        ACT_GATE_X0_F = float(args.act_gate_x0_f)
    if args.act_gate_k_f is not None:
        ACT_GATE_K_F = float(args.act_gate_k_f)
    if args.ia_k_force_f is not None:
        IA_K_FORCE_F = float(args.ia_k_force_f)
    if args.ia_k_stretch_f is not None:
        IA_K_STRETCH_F = float(args.ia_k_stretch_f)

    # inhibitory compensation: direct overrides first, then the gain multiplies them
    if args.w_ine2rgf is not None: W_INE2RGF = float(args.w_ine2rgf)
    if args.w_inf2rge is not None: W_INF2RGE = float(args.w_inf2rge)
    if args.i_e_rgf   is not None: I_E_RGF   = float(args.i_e_rgf)
    if args.w_rg_rec_f is not None: W_RG_REC_F = float(args.w_rg_rec_f)
    _g  = float(args.inh_comp)
    _gf = float(args.inh_comp_f)
    W_RG2INE  *= _g;   W_INE2RGF *= _g
    W_RG2INF  *= _gf;  W_INF2RGE *= _gf

    # ablations
    ablation_tag=[]
    if args.ablate_ia_loop: W_IA2IN=0.0; ablation_tag.append("noIa")
    if args.ablate_asym:    W_INF2RGE=-28.0; W_INE2RGF=-28.0; ablation_tag.append("symInh")
    if args.ablate_comm:    W_COMM_F_INH=0.0; W_COMM_E_INH=0.0; ablation_tag.append("noComm")
    BS_RATE_BASE_HZ=float(args.bs_base_hz)
    ENFORCE_TONIC_BS=bool(args.enforce_tonic_bs)
    IA_FEEDBACK_GAIN=float(args.ia_feedback_gain)
    CUT_FEEDBACK_GAIN=float(args.cut_feedback_gain)
    PACED_GAIT=bool(getattr(args,"paced_gait",False))
    LAMBDA=0.001
    if args.stdp_lambda: LAMBDA=float(args.stdp_lambda)

    # sweep
    def _pp(s):
        out=[]
        for item in (s or "").strip().split(","):
            item=item.strip()
            if not item: continue
            mu_s,cv_s=item.split(":",1); out.append((float(mu_s),float(cv_s)))
        return out
    sweep_pairs=_pp(getattr(args,"sweep_pairs",""))
    sweep_active=len(sweep_pairs)>0
    sweep_idx=int(getattr(args,"sweep_run_idx",-1))
    if sweep_active and sweep_idx<0:
        sweep_idx=int(os.environ.get("SLURM_ARRAY_TASK_ID","0"))
    sweep_mu=sweep_cv=sweep_sigma=None
    run_seed=int(getattr(args,"seed",12345))
    if sweep_active:
        sweep_mu,sweep_cv=sweep_pairs[sweep_idx]
        sweep_sigma=float(sweep_cv)*float(sweep_mu)
        run_seed=int(getattr(args,"seed",12345))+int(sweep_idx)*10007
        np.random.seed(run_seed)
        outdir=str(getattr(args,"outdir","."))
        tag=str(getattr(args,"tag","cpg_fixedw"))
        rn=str(getattr(args,"run_name","")).strip()
        base=rn if rn else f"cpg_{tag}_idx{sweep_idx:02d}_mu{sweep_mu:05.2f}_cv{sweep_cv:05.2f}_seed{run_seed}"
        if not base.endswith(".h5"): base+=".h5"
        if str(getattr(args,"out","")).strip() in ("","cpg_run.h5"):
            args.out=os.path.join(outdir,base)

    try: nest.set_verbosity(str(args.nest_verbosity))
    except: pass
    if args.long_run:
        if args.simulate_chunk_ms<100.: args.simulate_chunk_ms=100.  # do NOT lower: rate feeds the gate, not just the plot
        if args.rate_update_ms<100.:    args.rate_update_ms=100.      # 20ms bins -> spike-count noise -> corr(RG) -0.91 to -0.55
        if args.weight_sample_ms<1000.: args.weight_sample_ms=1000.
        if args.print_every<200:        args.print_every=200

    SIM_MS=float(args.sim_ms); DT_MS=float(args.dt_ms)
    PHASE_MS=SIM_MS/int(N_PHASES); CHUNK_MS=float(args.simulate_chunk_ms)
    RES_MS=float(args.resolution_ms)
    def q(x):
        if RES_MS<=0.: return float(x)
        return int(round(float(x)/RES_MS))*RES_MS
    CHUNK_MS=q(CHUNK_MS)
    if CHUNK_MS<=0.: raise ValueError("chunk<=0")
    if CHUNK_MS<DT_MS: raise ValueError(f"chunk {CHUNK_MS}<dt {DT_MS}")
    rate_every=max(1,int(round(float(args.rate_update_ms)/CHUNK_MS)))

    if PACED_GAIT:
        SP=q(float(args.step_period_ms)); STANCE_FRAC=float(args.stance_fraction)
        HALF=q(SP/2.); STANCE=q(SP*STANCE_FRAC); SWING=q(HALF-STANCE)
        NIG=int(args.n_ia_groups)
        IEH=list(args.ia_ext_hz)
        while len(IEH)<NIG: IEH.append(float(IEH[-1]))
        IEF=float(args.ia_ext_f_hz); SUB=q(STANCE/max(1,NIG))
        NHC=max(2,int(np.ceil(SIM_MS/HALF)))
        # N confirmed the baseline itself uses these under paced gait
        # (verified against the real source model) -- keeping them here.
        TAU_ACT_RISE_MS=TAU_ACT_DECAY_MS=40.
        TAU_FORCE_RISE_MS=TAU_FORCE_DECAY_MS=80.

    nest.ResetKernel()
    # seed the NEST kernel too, not just numpy, or --seed does not control connectivity
    nest.SetKernelStatus({"resolution":RES_MS,"local_num_threads":int(args.threads),
                          "print_time":False,"rng_seed":int(run_seed)})
    print(f"[varAF] NEST kernel rng_seed={run_seed} threads={args.threads}")
    dm=str(getattr(args,"delay_model","fixed")); sp=str(getattr(args,"species","rat"))
    dj=float(getattr(args,"delay_jitter_ms",0.)); ds=float(getattr(args,"delay_scale",1.))
    DL={k:_mk_delay(dm,sp,k,fb,RES_MS,dj,ds) for k,fb in [
        ("cut_to_rg",DELAY_MS),("bs_to_rg",DELAY_MS),("base_to_rg",DELAY_MS),
        ("rg_to_m",DELAY_MS),("m_to_mus",DELAY_MS),("ia_path",DELAY_MS),
        ("rg_rec",DELAY_MS),("rg_recip",DELAY_RECIP_MS),
        ("motor_e2f",DELAY_MOTOR_RECIP_E2F_MS),("motor_f2e",DELAY_MOTOR_RECIP_F2E_MS),
        ("commissural",DELAY_COMM_MS)]}

    try: os.makedirs(os.path.dirname(str(args.out)),exist_ok=True)
    except: pass
    if "SLURM_PROCID" in os.environ:
        rank=int(os.environ.get("SLURM_PROCID","0")); nproc=int(os.environ.get("SLURM_NTASKS","1"))
    else:
        rank=getattr(nest,"Rank",lambda:0)(); nproc=getattr(nest,"NumProcesses",lambda:1)()

    if rank==0:
        print("="*65)
        print(f"FIXED-WEIGHTS VARIANT Q  k={_k}  m={_m}  seed={run_seed}")
        print(f"  STDP DISABLED: CUT/BS->RG all static_synapse weight={W0_IN}")
        print(f"  RG circuit = N's validated state: I_E_RGF={I_E_RGF} W_INE2RGF={W_INE2RGF}")
        print(f"  [varX] inh_comp={args.inh_comp} inh_comp_f={args.inh_comp_f} -> W_RG2INE={W_RG2INE} W_INE2RGF={W_INE2RGF} W_RG2INF={W_RG2INF} W_INF2RGE={W_INF2RGE}")
        print(f"  W_RG_REC_F={W_RG_REC_F}  P_RECIP_E={P_RG_RECIP_E}  P_RECIP_F={P_RG_RECIP_F}")
        print(f"  *** NEW: split activation-gate threshold ***")
        print(f"  ACT_GATE_X0_E={ACT_GATE_X0_E} (K={ACT_GATE_K})  ACT_GATE_X0_F={ACT_GATE_X0_F} (K={ACT_GATE_K_F})  [gate-recalibrated to measured RG-F range]")
        print(f"  IA_K_FORCE_F={IA_K_FORCE_F}  IA_K_STRETCH_F={IA_K_STRETCH_F}  [Ia-F feedback RESTORED, was 0.0/0.0]")
        print(f"  TAU_ACT={TAU_ACT_RISE_MS}/{TAU_ACT_DECAY_MS}  TAU_FORCE={TAU_FORCE_RISE_MS}/{TAU_FORCE_DECAY_MS}")
        print("="*65)

    # ── Build network ───────────────────────────────────────
    leg={}
    for side in LEGS:
        cut_pg=nest.Create("poisson_generator",N_CUT)
        cut_in=nest.Create("parrot_neuron",N_CUT)
        nest.Connect(cut_pg,cut_in,conn_spec={"rule":"one_to_one"})
        nest.SetStatus(cut_pg,{"rate":CUT_RATE_OFF_HZ})

        bs_pg_e=nest.Create("spike_generator",N_BS)
        bs_in_e=nest.Create("parrot_neuron",N_BS)
        nest.Connect(bs_pg_e,bs_in_e,conn_spec={"rule":"one_to_one"})
        bs_pg_f=nest.Create("spike_generator",N_BS)
        bs_in_f=nest.Create("parrot_neuron",N_BS)
        nest.Connect(bs_pg_f,bs_in_f,conn_spec={"rule":"one_to_one"})

        period=1000./float(BS_REGULAR_HZ)
        bt=np.arange(RES_MS,SIM_MS+1e-9,period)
        offs=np.linspace(0.,period,int(N_BS),endpoint=False)
        te=[]; tf=[]
        for off in offs:
            ta=bt+float(off)
            if BS_REGULAR_JITTER_MS>0.:
                j=float(BS_REGULAR_JITTER_MS)
                ta=np.clip(ta+np.random.uniform(-j,j,size=ta.shape),RES_MS,SIM_MS)
                ta=np.round(ta/RES_MS)*RES_MS; ta=ta[ta>=RES_MS]; ta=np.unique(np.sort(ta))
            te.append(ta.tolist()); tf.append(ta.tolist())
        nest.SetStatus(bs_pg_e,[{"spike_times":st} for st in te])
        nest.SetStatus(bs_pg_f,[{"spike_times":st} for st in tf])

        base_pg=nest.Create("poisson_generator",N_BS); base_in=nest.Create("parrot_neuron",N_BS)
        nest.Connect(base_pg,base_in,conn_spec={"rule":"one_to_one"}); nest.SetStatus(base_pg,{"rate":2.0})

        ia_pg_e=nest.Create("poisson_generator",N_IA_E); ia_in_e=nest.Create("parrot_neuron",N_IA_E)
        nest.Connect(ia_pg_e,ia_in_e,conn_spec={"rule":"one_to_one"}); nest.SetStatus(ia_pg_e,{"rate":IA_BASE_HZ_E})
        ia_pg_f=nest.Create("poisson_generator",N_IA_F); ia_in_f=nest.Create("parrot_neuron",N_IA_F)
        nest.Connect(ia_pg_f,ia_in_f,conn_spec={"rule":"one_to_one"}); nest.SetStatus(ia_pg_f,{"rate":IA_BASE_HZ_F})

        ia_ext_pg_e=[]; ia_ext_pg_f=None
        if PACED_GAIT:
            np_=max(1,N_IA_E//NIG)
            for _ in range(NIG):
                pg=nest.Create("poisson_generator",np_); pn=nest.Create("parrot_neuron",np_)
                nest.Connect(pg,pn,conn_spec={"rule":"one_to_one"}); nest.SetStatus(pg,{"rate":0.})
                ia_ext_pg_e.append(pg)
            if IEF>0.:
                ia_ext_pg_f=nest.Create("poisson_generator",N_IA_F); nest.SetStatus(ia_ext_pg_f,{"rate":0.})

        rg_e=nest.Create("izhikevich",N_RG_E); rg_f=nest.Create("izhikevich",N_RG_F)
        m_e =nest.Create("izhikevich",N_MOTOR_E); m_f=nest.Create("izhikevich",N_MOTOR_F)
        rec_rge=nest.Create("spike_recorder"); nest.Connect(rg_e,rec_rge)
        rec_rgf=nest.Create("spike_recorder"); nest.Connect(rg_f,rec_rgf)
        ia_int_e=nest.Create("izhikevich",N_IA_INT); ia_int_f=nest.Create("izhikevich",N_IA_INT)
        in_e=nest.Create("izhikevich",N_INE); in_f=nest.Create("izhikevich",N_INF)
        rec_ine=nest.Create("spike_recorder"); nest.Connect(in_e,rec_ine)
        rec_inf=nest.Create("spike_recorder"); nest.Connect(in_f,rec_inf)
        rec_iainte=nest.Create("spike_recorder"); nest.Connect(ia_int_e,rec_iainte)
        rec_iaintf=nest.Create("spike_recorder"); nest.Connect(ia_int_f,rec_iaintf)
        for pop in (rg_e,rg_f,m_e,m_f): nest.SetStatus(pop,izh_params)
        for pop in (ia_int_e,ia_int_f,in_e,in_f): nest.SetStatus(pop,izh_inh_params)
        bias=+LEFT_RIGHT_BIAS_IE if side=="L" else -LEFT_RIGHT_BIAS_IE
        nest.SetStatus(rg_e,{"V_m":-65.,"U_m":0.2*(-65.),"I_e":I_E_RGE})
        nest.SetStatus(rg_f,{"a":RGF_A,"b":RGF_B,"c":RGF_C,"d":RGF_D,"V_m":-65.,"U_m":RGF_B*(-65.),"I_e":I_E_RGF+bias})
        nest.SetStatus(m_e,{"V_m":-65.,"U_m":0.2*(-65.),"I_e":I_E_MOTOR})
        nest.SetStatus(m_f,{"V_m":-65.,"U_m":0.2*(-65.),"I_e":I_E_MOTOR})
        mus_e=nest.Create("parrot_neuron",N_MUS_E); mus_f=nest.Create("parrot_neuron",N_MUS_F)
        rec_muse=nest.Create("spike_recorder"); nest.Connect(mus_e,rec_muse)
        rec_musf=nest.Create("spike_recorder"); nest.Connect(mus_f,rec_musf)
        leg[side]=dict(cut_pg=cut_pg,cut_in=cut_in,bs_pg_e=bs_pg_e,bs_in_e=bs_in_e,
            bs_pg_f=bs_pg_f,bs_in_f=bs_in_f,base_pg=base_pg,base_in=base_in,
            ia_pg_e=ia_pg_e,ia_in_e=ia_in_e,ia_pg_f=ia_pg_f,ia_in_f=ia_in_f,
            ia_ext_pg_e=ia_ext_pg_e,ia_ext_pg_f=ia_ext_pg_f,
            rg_e=rg_e,rg_f=rg_f,m_e=m_e,m_f=m_f,
            ia_int_e=ia_int_e,ia_int_f=ia_int_f,in_e=in_e,in_f=in_f,
            mus_e=mus_e,mus_f=mus_f,rec_muse=rec_muse,rec_musf=rec_musf,
            rec_rge=rec_rge,rec_rgf=rec_rgf,rec_ine=rec_ine,rec_inf=rec_inf,
            rec_iainte=rec_iainte,rec_iaintf=rec_iaintf)

    # ── plastic input synapse model (STDP off when lambda==0) ──
    _stdp_lambda = float(getattr(args,"stdp_lambda",0.0))
    _use_stdp = _stdp_lambda > 0.0
    if _use_stdp:
        nest.SetDefaults("stdp_synapse", {
            "lambda": _stdp_lambda, "alpha": 1.0,
            "mu_plus": 1.0, "mu_minus": 1.0,
            "Wmax": float(getattr(args,"stdp_wmax",40.0)), "tau_plus": 20.0})
        print(f"  [stdp] plastic cut->RG-E, bs->RG-E, bs->RG-F  lambda={_stdp_lambda} "
              f"Wmax={getattr(args,'stdp_wmax',40.0)}")

    # ── Connect  (reciprocal inhibition ALWAYS static; inputs optionally plastic) ──
    def S(src,tgt,p,w,d): nest.Connect(src,tgt,conn_spec={"rule":"pairwise_bernoulli","p":p},
                                        syn_spec={"synapse_model":"static_synapse","weight":w,"delay":d})
    def SIN(src,tgt,p,w,d):
        # sensory/descending input to RG: plastic if STDP on, else static
        model = "stdp_synapse" if _use_stdp else "static_synapse"
        nest.Connect(src,tgt,conn_spec={"rule":"pairwise_bernoulli","p":p},
                     syn_spec={"synapse_model":model,"weight":w,"delay":d})
    for side in LEGS:
        L=leg[side]
        SIN(L["cut_in"],L["rg_e"],P_IN_STDP,W0_IN,DL["cut_to_rg"])
        SIN(L["bs_in_e"],L["rg_e"],P_IN_STDP,W0_IN,DL["bs_to_rg"])
        SIN(L["bs_in_f"],L["rg_f"],P_IN_STDP,W0_IN,DL["bs_to_rg"])
        S(L["base_in"],L["rg_e"],BASE_DRIVE_P,BASE_DRIVE_W,DL["base_to_rg"])
        S(L["base_in"],L["rg_f"],BASE_DRIVE_P,BASE_DRIVE_W,DL["base_to_rg"])
        S(L["rg_e"],L["m_e"],P_IN_STDP,W0_RM,DL["rg_to_m"])
        S(L["rg_f"],L["m_f"],P_IN_STDP,W0_RM,DL["rg_to_m"])
        S(L["m_e"],L["m_f"],P_MOTOR_RECIP,W_MOTOR_RECIP,DL["motor_e2f"])
        S(L["m_f"],L["m_e"],P_MOTOR_RECIP,W_MOTOR_RECIP,DL["motor_f2e"])
        S(L["m_e"],L["mus_e"],P_M2MUS,W_M2MUS,DL["m_to_mus"])
        S(L["m_f"],L["mus_f"],P_M2MUS,W_M2MUS,DL["m_to_mus"])
        S(L["ia_in_e"],L["ia_int_e"],IA2RG_P,W_IA_IN2INT,DL["ia_path"])
        S(L["ia_int_e"],L["m_f"],    IA2RG_P,W_IA_INT2ANT,DL["ia_path"])
        S(L["ia_in_f"],L["ia_int_f"],IA2RG_P,W_IA_IN2INT,DL["ia_path"])
        S(L["ia_int_f"],L["m_e"],    IA2RG_P,W_IA_INT2ANT,DL["ia_path"])
        S(L["ia_in_e"],L["in_e"],P_IA2IN,W_IA2IN,DL["ia_path"])
        S(L["ia_in_f"],L["in_f"],P_IA2IN,W_IA2IN,DL["ia_path"])
        if PACED_GAIT:
            for pn in L["ia_ext_pg_e"]: S(pn,L["in_e"],P_IA2IN,W_IA2IN,DL["ia_path"])
            if L["ia_ext_pg_f"] is not None:
                S(L["ia_ext_pg_f"],L["rg_f"],P_IA2IN,W_FLEX_AFF2RGF,DL["ia_path"])
                S(L["ia_ext_pg_f"],L["in_f"],P_IA2IN,W_IA2IN,DL["ia_path"])
        S(L["rg_e"],L["rg_e"],P_RG_REC,W_RG_REC_E,DL["rg_rec"])
        S(L["rg_f"],L["rg_f"],P_RG_REC,W_RG_REC_F,DL["rg_rec"])
        S(L["rg_f"],L["in_f"],P_RG_RECIP_F,W_RG2INF,DL["rg_recip"])
        S(L["in_f"],L["rg_e"],P_RG_RECIP_F,W_INF2RGE,DL["rg_recip"])
        S(L["rg_e"],L["in_e"],P_RG_RECIP_E,W_RG2INE,DL["rg_recip"])
        S(L["in_e"],L["rg_f"],P_RG_RECIP_E,W_INE2RGF,DL["rg_recip"])
        S(L["cut_in"],L["in_e"],P_CUT2INE,W_CUT2INE,DL["cut_to_rg"])

    if ENABLE_COMMISSURAL:
        LL=leg["L"]; RR=leg["R"]
        for src,tgt,p,w in [(LL["rg_f"],RR["rg_f"],P_COMM_F,W_COMM_F_INH),
                             (RR["rg_f"],LL["rg_f"],P_COMM_F,W_COMM_F_INH),
                             (LL["rg_e"],RR["rg_e"],P_COMM_E,W_COMM_E_INH),
                             (RR["rg_e"],LL["rg_e"],P_COMM_E,W_COMM_E_INH)]:
            S(src,tgt,p,w,DL["commissural"])

    # Bio-plausible weight heterogeneity
    scv=float(getattr(args,"static_weight_cv",0.0) or 0.0)
    if scv>0.:
        try:
            sc=nest.GetConnections(synapse_model="static_synapse")
            if sc is not None and len(sc)>0:
                w=np.asarray(nest.GetStatus(sc,"weight"),dtype=float)
                sig=float(np.sqrt(np.log(1.+scv*scv))); mu=-0.5*sig*sig
                fac=np.random.default_rng(int(args.seed)+777).lognormal(mu,sig,size=w.size)
                nest.SetStatus(sc,[{"weight":float(wi*fi)} for wi,fi in zip(w,fac)])
        except Exception as e:
            if rank==0: print(f"[BioPlaus] skipped: {e}")

    stats_syn={"static_total":safe_len("static_synapse")}
    if rank==0:
        print(f"[Stats] static synapses={stats_syn['static_total']}  NO STDP synapses.")

    plastic_keys=[]

    times=[]; logs={s:dict(bs_e=[],bs_f=[],mus_e=[],mus_f=[],
        rge=[],rgf=[],ine=[],inf=[],iaint_e=[],iaint_f=[],
        act_e=[],act_f=[],force_e=[],force_f=[],len_e=[],len_f=[],
        ia_e=[],ia_f=[]) for s in LEGS}
    state={s:dict(act_e=0.,act_f=0.,force_e=0.,force_f=0.,len_e=L0,len_f=L0,
                  last_muse=0,last_musf=0,last_rge=0,last_rgf=0,
                  last_ine=0,last_inf=0,last_iainte=0,last_iaintf=0) for s in LEGS}

    def update_leg(side,t_ms,dt_ms_,cut_frac,do_rate):
        dt_s=float(dt_ms_)/1000.; L=leg[side]; St=state[side]; P=logs[side]
        r_e=r_f=clamp(BS_REGULAR_HZ,0.,BS_REGULAR_HZ)
        P["bs_e"].append(r_e); P["bs_f"].append(r_f)
        sp_e,ce=new_spikes(L["rec_muse"],St["last_muse"]); St["last_muse"]=ce
        sp_f,cf=new_spikes(L["rec_musf"],St["last_musf"]); St["last_musf"]=cf
        sp_re,cre=new_spikes(L["rec_rge"],St["last_rge"]); St["last_rge"]=cre
        sp_rf,crf=new_spikes(L["rec_rgf"],St["last_rgf"]); St["last_rgf"]=crf
        dts=max(1e-9,dt_s)
        fi_e=max(1.,float(N_MOTOR_E)*float(P_M2MUS)); fi_f=max(1.,float(N_MOTOR_F)*float(P_M2MUS))
        r_me=((sp_e/max(1,N_MUS_E))/dts)/fi_e; r_mf=((sp_f/max(1,N_MUS_F))/dts)/fi_f
        P["mus_e"].append(r_me); P["mus_f"].append(r_mf)
        r_rge=(sp_re/max(1,N_RG_E))/dts; r_rgf=(sp_rf/max(1,N_RG_F))/dts
        P["rge"].append(r_rge); P["rgf"].append(r_rgf)
        sp_ie,cie=new_spikes(L["rec_ine"],St["last_ine"]); St["last_ine"]=cie
        sp_if,cif=new_spikes(L["rec_inf"],St["last_inf"]); St["last_inf"]=cif
        sp_ae,cae=new_spikes(L["rec_iainte"],St["last_iainte"]); St["last_iainte"]=cae
        sp_af,caf=new_spikes(L["rec_iaintf"],St["last_iaintf"]); St["last_iaintf"]=caf
        P["ine"].append((sp_ie/max(1,N_INE))/dts); P["inf"].append((sp_if/max(1,N_INF))/dts)
        P["iaint_e"].append((sp_ae/max(1,N_IA_INT))/dts); P["iaint_f"].append((sp_af/max(1,N_IA_INT))/dts)
        ref=100.
        # VARIANT Q: independent gate thresholds for E and F
        de=1./(1.+np.exp(-ACT_GATE_K*(float(r_rge)/ref-ACT_GATE_X0_E)))
        df=1./(1.+np.exp(-ACT_GATE_K_F*(float(r_rgf)/ref-ACT_GATE_X0_F)))
        _GATE_DIAG["de"].append(float(de)); _GATE_DIAG["df"].append(float(df))
        ae=ACT_MAX*(1.-np.exp(-ACT_SAT_K*float(r_me)))*de
        af=ACT_MAX*(1.-np.exp(-ACT_SAT_K*float(r_mf)))*df
        tr=TAU_ACT_RISE_MS/1000.; td=TAU_ACT_DECAY_MS/1000.
        ke=1.-np.exp(-dts/max(1e-9,tr if ae>St["act_e"] else td))
        kf=1.-np.exp(-dts/max(1e-9,tr if af>St["act_f"] else td))
        St["act_e"]+=ke*(ae-St["act_e"]); St["act_f"]+=kf*(af-St["act_f"])
        St["act_e"]=clamp(St["act_e"],0.,ACT_MAX); St["act_f"]=clamp(St["act_f"],0.,ACT_MAX)
        tfe=FORCE_MAX*(1.-np.exp(-FORCE_SAT_K*St["act_e"]))
        tff=FORCE_MAX*(1.-np.exp(-FORCE_SAT_K*St["act_f"]))
        tr2=TAU_FORCE_RISE_MS/1000.; td2=TAU_FORCE_DECAY_MS/1000.
        kfe=1.-np.exp(-dts/max(1e-9,tr2 if tfe>St["force_e"] else td2))
        kff=1.-np.exp(-dts/max(1e-9,tr2 if tff>St["force_f"] else td2))
        St["force_e"]+=kfe*(tfe-St["force_e"]); St["force_f"]+=kff*(tff-St["force_f"])
        St["force_e"]=clamp(St["force_e"],0.,FORCE_MAX); St["force_f"]=clamp(St["force_f"],0.,FORCE_MAX)
        tl=TAU_LENGTH_MS/1000.; kl=1.-np.exp(-dts/max(1e-9,tl))
        St["len_e"]+=kl*(L0-St["len_e"]); St["len_f"]+=kl*(L0-St["len_f"])
        St["len_e"]-=SHORTEN_GAIN*St["force_e"]*dt_s; St["len_f"]-=SHORTEN_GAIN*St["force_f"]*dt_s
        if cut_frac>0.: St["len_e"]+=STRETCH_GAIN*cut_frac*dt_s
        St["len_e"]=clamp(St["len_e"],L_MIN,L_MAX); St["len_f"]=clamp(St["len_f"],L_MIN,L_MAX)
        se=max(0.,St["len_e"]-L0); sf=max(0.,St["len_f"]-L0)
        ia_e=IA_FEEDBACK_GAIN*clamp(IA_BASE_HZ_E+IA_K_FORCE_E*St["force_e"]+IA_K_STRETCH_E*se,0.,IA_RATE_MAX_HZ)
        ia_f=IA_FEEDBACK_GAIN*clamp(IA_BASE_HZ_F+IA_K_FORCE_F*St["force_f"]+IA_K_STRETCH_F*sf,0.,IA_RATE_MAX_HZ)
        if do_rate: nest.SetStatus(L["ia_pg_e"],{"rate":ia_e}); nest.SetStatus(L["ia_pg_f"],{"rate":ia_f})
        P["act_e"].append(St["act_e"]); P["act_f"].append(St["act_f"])
        P["force_e"].append(St["force_e"]); P["force_f"].append(St["force_f"])
        P["len_e"].append(St["len_e"]); P["len_f"].append(St["len_f"])
        P["ia_e"].append(ia_e); P["ia_f"].append(ia_f)

    done=0; t_ms=0.; t0=time.time(); sa=0.; ba=0.
    def run_win(win,caf):
        nonlocal done,t_ms,sa,ba
        nc=int(win//CHUNK_MS); tail=q(win-nc*CHUNK_MS); nct=nc+(1 if tail>1e-9 else 0)
        for ci in range(nct):
            cur=q(CHUNK_MS if ci<nc else tail)
            if cur<=0.: continue
            ts=time.perf_counter(); nest.Simulate(cur); sa+=time.perf_counter()-ts
            t_ms+=cur; done+=1; tb=time.perf_counter()
            dr=(done%rate_every==0)
            for s in LEGS: update_leg(s,t_ms,cur,caf,dr)
            times.append(t_ms); ba+=time.perf_counter()-tb

    if PACED_GAIT:
        for hc in range(NHC):
            if t_ms>=SIM_MS: break
            stance="L" if hc%2==0 else "R"; swing="R" if hc%2==0 else "L"
            nest.SetStatus(leg[swing]["cut_pg"],{"rate":CUT_RATE_OFF_HZ})
            for g in leg[swing]["ia_ext_pg_e"]: nest.SetStatus(g,{"rate":0.})
            if leg[swing]["ia_ext_pg_f"] is not None: nest.SetStatus(leg[swing]["ia_ext_pg_f"],{"rate":IEF})
            nest.SetStatus(leg[stance]["cut_pg"],{"rate":CUT_FEEDBACK_GAIN*CUT_RATE_ON_HZ})
            if leg[stance]["ia_ext_pg_f"] is not None: nest.SetStatus(leg[stance]["ia_ext_pg_f"],{"rate":0.})
            for gi in range(NIG):
                nest.SetStatus(leg[stance]["ia_ext_pg_e"][gi],{"rate":IEH[gi]})
                if gi>0: nest.SetStatus(leg[stance]["ia_ext_pg_e"][gi-1],{"rate":0.})
                sr=min(SUB,SIM_MS-t_ms)
                if sr<=0.: break
                run_win(sr,1.)
            for g in leg[stance]["ia_ext_pg_e"]: nest.SetStatus(g,{"rate":0.})
            if SWING>0.:
                nest.SetStatus(leg[stance]["cut_pg"],{"rate":CUT_RATE_OFF_HZ})
                swr=min(SWING,SIM_MS-t_ms)
                if swr>0.: run_win(swr,0.)
            if rank==0 and hc%max(1,args.print_every)==0:
                print(f"[Paced] hc {hc+1}/{NHC} stance={stance} t={t_ms:.0f}/{SIM_MS:.0f}")
    else:
        cc=max(1,int(N_CUT/N_PHASES))
        for phase in range(N_PHASES):
            for s in LEGS: nest.SetStatus(leg[s]["cut_pg"],{"rate":CUT_RATE_OFF_HZ})
            st=phase*cc; en=min(N_CUT,(phase+1)*cc)
            for s in LEGS: nest.SetStatus(leg[s]["cut_pg"][st:en],{"rate":CUT_FEEDBACK_GAIN*CUT_RATE_ON_HZ})
            caf=float(en-st)/float(N_CUT)
            nc=int(PHASE_MS//CHUNK_MS); tail=q(PHASE_MS-nc*CHUNK_MS); nct=nc+(1 if tail>1e-9 else 0)
            for lc in range(nct):
                cur=q(CHUNK_MS if lc<nc else tail)
                if cur<=0.: continue
                ts=time.perf_counter(); nest.Simulate(cur); sa+=time.perf_counter()-ts
                t_ms+=cur; done+=1; tb=time.perf_counter()
                dr=(done%rate_every==0)
                for s in LEGS: update_leg(s,t_ms,cur,caf,dr)
                times.append(t_ms); ba+=time.perf_counter()-tb
                if rank==0 and done%int(args.print_every)==0:
                    print(f"[Sim] phase {phase+1}/{N_PHASES} chunk {done} t={t_ms:.0f}")

    if rank==0: print(f"[Done] {time.time()-t0:.1f}s  out={args.out}")
    if rank!=0: return

    # ── Write HDF5 ──────────────────────────────────────────
    # attrs wrapped in try/except so a metadata failure cannot abort the dataset write
    with h5py.File(args.out,"w") as h5:
      try:
        h5.attrs["created_utc"]=datetime.utcnow().isoformat()+"Z"
        h5.attrs["nest_version"]=str(nest.__version__)
        h5.attrs["sim_ms"]=SIM_MS; h5.attrs["dt_ms"]=CHUNK_MS
        h5.attrs["resolution_ms"]=float(args.resolution_ms)
        h5.attrs["paced_gait"]=bool(PACED_GAIT)
        h5.attrs["ablation_tag"]=",".join(ablation_tag) if ablation_tag else "baseline"
        h5.attrs["delay_model"]=str(getattr(args,"delay_model","fixed"))
        h5.attrs["species"]=str(getattr(args,"species","rat"))
        h5.attrs["stdp_lambda"]=float(getattr(args,"stdp_lambda",0.0))
        h5.attrs["fixed_weights"]=(float(getattr(args,"stdp_lambda",0.0))==0.0)
        h5.attrs["variant"]="T_k025_m040_fixedweights_gate_recalibrated_IaF_restored"
        h5.attrs["k_connectivity"]=float(_k)
        h5.attrs["m_neurons"]=float(_m)
        h5.attrs["I_E_RGF"]=float(I_E_RGF)
        h5.attrs["W_RG_REC_F"]=float(W_RG_REC_F)
        h5.attrs["W_INE2RGF"]=float(W_INE2RGF)
        h5.attrs["W_INF2RGE"]=float(W_INF2RGE)
        h5.attrs["W_RG2INE"]=float(W_RG2INE)
        h5.attrs["W_RG2INF"]=float(W_RG2INF)
        h5.attrs["ACT_GATE_X0_E_used"]=float(ACT_GATE_X0_E)
        h5.attrs["ACT_GATE_K_E_used"]=float(ACT_GATE_K)
        h5.attrs["ACT_GATE_X0_F_used"]=float(ACT_GATE_X0_F)
        h5.attrs["ACT_GATE_K_F_used"]=float(ACT_GATE_K_F)
        # note: trace dict is logs[side] here; P only exists inside update_leg()
        _dg = logs["L"]
        for _nm in ("rge","rgf"):
            _v=np.asarray(_dg.get(_nm,[]),dtype=float)
            if _v.size:
                _p1,_p99=float(np.percentile(_v,1)),float(np.percentile(_v,99))
                h5.attrs[f"RG_{_nm.upper()}_min"]=_p1
                h5.attrs[f"RG_{_nm.upper()}_max"]=_p99
                h5.attrs[f"RG_{_nm.upper()}_ptp"]=_p99-_p1
                print(f"  [varY] RG-{_nm[-1].upper()} rate Hz: p1={_p1:.1f} p99={_p99:.1f} ptp={_p99-_p1:.1f}")
        for _nm,_key in (("de","act_e"),("df","act_f")):
            _v=np.asarray(_GATE_DIAG[_nm],dtype=float)
            if _v.size:
                _p1,_p99=float(np.percentile(_v,1)),float(np.percentile(_v,99))
                h5.attrs[f"GATE_{_nm}_min"]=_p1
                h5.attrs[f"GATE_{_nm}_max"]=_p99
                h5.attrs[f"GATE_{_nm}_ptp"]=_p99-_p1
                print(f"  [varY] gate {_nm}: p1={_p1:.3f} p99={_p99:.3f} ptp={_p99-_p1:.3f}")
        for _nm in ("force_e","force_f"):
            _v=np.asarray(_dg.get(_nm,[]),dtype=float)
            if _v.size:
                print(f"  [varY] {_nm}: p1={np.percentile(_v,1):.2f} p99={np.percentile(_v,99):.2f}")
        h5.attrs["INH_COMP"]=float(args.inh_comp)
        h5.attrs["INH_COMP_F"]=float(args.inh_comp_f)
        h5.attrs["k_conn"]=float(_k)
        h5.attrs["m_neuron"]=float(_m)
        h5.attrs["N_RG_E"]=int(N_RG_E); h5.attrs["N_INE"]=int(N_INE)
        h5.attrs["P_RG_RECIP_E"]=float(P_RG_RECIP_E)
        h5.attrs["P_RG_RECIP_F"]=float(P_RG_RECIP_F)
        h5.attrs["ACT_GATE_X0_E"]=float(ACT_GATE_X0_E)
        h5.attrs["ACT_GATE_X0_F"]=float(ACT_GATE_X0_F)
        h5.attrs["ACT_GATE_K_E"]=float(ACT_GATE_K)
        h5.attrs["ACT_GATE_K_F"]=float(ACT_GATE_K_F)
        h5.attrs["IA_K_FORCE_F"]=float(IA_K_FORCE_F)
        h5.attrs["IA_K_STRETCH_F"]=float(IA_K_STRETCH_F)
        h5.attrs["TAU_ACT_RISE_MS"]=float(TAU_ACT_RISE_MS)
        h5.attrs["TAU_ACT_DECAY_MS"]=float(TAU_ACT_DECAY_MS)
        h5.attrs["TAU_FORCE_RISE_MS"]=float(TAU_FORCE_RISE_MS)
        h5.attrs["TAU_FORCE_DECAY_MS"]=float(TAU_FORCE_DECAY_MS)
        h5.attrs["freeze_bs_rg"]=True
        h5.attrs["stdp_ia_rg"]=False
        h5.attrs["sweep_active"]=bool(sweep_active)
        if sweep_active:
            h5.attrs["sweep_idx"]=int(sweep_idx)
            h5.attrs["winit_mu"]=float(sweep_mu)
            h5.attrs["winit_cv"]=float(sweep_cv)
            h5.attrs["seed"]=int(run_seed)
      except Exception as _e:
        print(f"[WARN] attribute/diagnostic block failed ({_e!r}) -- datasets still written")
      gs=h5.create_group("stats"); gs.attrs["static_total"]=int(stats_syn["static_total"])
      h5.create_dataset("times_ms",data=np.asarray(times,dtype=np.float32),compression="gzip")
      for side in LEGS:
            g=h5.create_group(f"leg_{side}")
            for key,arr in logs[side].items():
                g.create_dataset(key,data=np.asarray(arr,dtype=np.float32),compression="gzip")
            g.create_group("weights")  # empty — no STDP
    print(f"[HDF5] saved → {args.out}")

if __name__=="__main__": main()
