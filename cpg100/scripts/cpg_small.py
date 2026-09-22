#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone

import numpy as np

try:
    import h5py
except ImportError:
    h5py = None
try:
    import nest
except ImportError:
    nest = None

SCRIPT_VERSION = "cpg_small 1.3"
LEGS = ("L", "R")

# ══════════════════════════════════════════════════════════════════════
# Sizes
# ══════════════════════════════════════════════════════════════════════
SIZE_KEYS = ("RGE", "RGF", "INE", "INF", "IAINT", "MNE", "MNF")
# control cell (k=0.20, m=0.40): the reference for --preserve-input
REF_SIZES = dict(RGE=40, RGF=40, INE=20, INF=20, IAINT=20, MNE=40, MNF=40)
RELAY_N_DEFAULT = 40


def izh_per_leg(s):
    return s["RGE"] + s["RGF"] + s["INE"] + s["INF"] + 2 * s["IAINT"] + s["MNE"] + s["MNF"]


def parse_sizes(txt):
    vals = [int(v) for v in str(txt).replace(" ", "").split(",") if v != ""]
    if len(vals) != len(SIZE_KEYS):
        raise SystemExit(f"--sizes needs {len(SIZE_KEYS)} integers "
                         f"({','.join(SIZE_KEYS)}), got {txt!r}")
    if min(vals) < 1:
        raise SystemExit("--sizes: every population needs at least 1 neuron")
    return dict(zip(SIZE_KEYS, vals))


# ══════════════════════════════════════════════════════════════════════
# Model constants (verbatim from cpg_optimization.py)
# ══════════════════════════════════════════════════════════════════════
W_IA_IN2INT = 6.0
W_IA_INT2ANT = -10.0
W_IA2IN = 6.0
W_FLEX_AFF2RGF = 6.0
W_RG2INE, W_RG2INF = 12.0, 18.0
W_INE2RGF = -8.0            # Zhang 6:1
W_INF2RGE = -48.0
W_RG_REC_E, W_RG_REC_F = 4.0, 5.5
W_MOTOR_RECIP = -22.0
W_M2MUS = 1.0
W_CUT2INE = 6.0
W_COMM_F_INH, W_COMM_E_INH = -20.0, -8.0
W0_IN = 22.0
W0_RM = 30.0
BASE_DRIVE_W = 1.0

# connection probabilities at k = 1 (multiplied by --k-conn where the model does)
P_IN_STDP = 0.5             # CUT/BS -> RG and RG -> motor: not k-scaled
P0_RG_REC = 0.12
P0_RG_RECIP_F = 0.30
P0_RG_RECIP_E = 0.15
P0_MOTOR_RECIP = 0.25
P0_M2MUS = 0.80
P0_IA2IN = 0.25
P0_COMM_F = 0.22
P0_COMM_E = 0.10
P_CUT2INE = 0.30
IA2RG_P = 0.4
BASE_DRIVE_P = 0.10

CUT_RATE_ON_HZ = 100.0
CUT_RATE_OFF_HZ = 0.0
BS_REGULAR_HZ = 60.0
BS_REGULAR_JITTER_MS = 0.3
BASE_DRIVE_HZ = 2.0
LEFT_RIGHT_BIAS_IE = 0.12

DELAY_MS = 1.0
DELAY_RECIP_MS = 1.0
DELAY_MOTOR_RECIP_E2F_MS = 1.5
DELAY_MOTOR_RECIP_F2E_MS = 1.0
DELAY_COMM_MS = 1.0
DELAY_PRESETS = {
    "rat": {
        "cut_to_rg":   {"syn_delay_ms": 1.0, "length_m": 0.005, "velocity_mps": 5.0},
        "bs_to_rg":    {"syn_delay_ms": 1.0, "length_m": 0.010, "velocity_mps": 5.0},
        "base_to_rg":  {"syn_delay_ms": 1.0, "length_m": 0.010, "velocity_mps": 5.0},
        "rg_to_m":     {"syn_delay_ms": 1.0, "length_m": 0.005, "velocity_mps": 5.0},
        "m_to_mus":    {"syn_delay_ms": 1.0, "length_m": 0.005, "velocity_mps": 5.0},
        "ia_path":     {"syn_delay_ms": 1.0, "length_m": 0.010, "velocity_mps": 10.0},
        "rg_rec":      {"syn_delay_ms": 0.8, "length_m": 0.002, "velocity_mps": 5.0},
        "rg_recip":    {"syn_delay_ms": 1.0, "length_m": 0.004, "velocity_mps": 5.0},
        "motor_e2f":   {"syn_delay_ms": 1.0, "length_m": 0.004, "velocity_mps": 5.0},
        "motor_f2e":   {"syn_delay_ms": 1.0, "length_m": 0.004, "velocity_mps": 5.0},
        "commissural": {"syn_delay_ms": 1.0, "length_m": 0.006, "velocity_mps": 5.0},
    },
}

izh_params = dict(a=0.02, b=0.2, c=-65.0, d=8.0, V_th=30.0, V_min=-120.0)
izh_inh_params = dict(a=0.1, b=0.2, c=-65.0, d=2.0, V_th=30.0, V_min=-120.0)
I_E_RGE, I_E_RGF, I_E_MOTOR = 1.0, 0.9, 1.0
RGF_A, RGF_B, RGF_C, RGF_D = 0.02, 0.2, -55.0, 4.0

# muscle / force / Ia model (unchanged; paced-gait filter values)
TAU_ACT_MS = 40.0
TAU_FORCE_MS = 80.0
ACT_MAX, ACT_SAT_K = 1.2, 0.02
FORCE_MAX, FORCE_SAT_K = 25.0, 1.0
TAU_LENGTH_MS, L0, L_MIN, L_MAX = 260.0, 1.0, 0.5, 2.0
SHORTEN_GAIN, STRETCH_GAIN = 0.010, 0.35
IA_BASE_HZ_E, IA_K_FORCE_E, IA_K_STRETCH_E = 10.0, 6.0, 250.0
IA_BASE_HZ_F, IA_K_FORCE_F, IA_K_STRETCH_F = 5.0, 6.0, 250.0
IA_RATE_MAX_HZ = 500.0
RATE_REF_HZ = 100.0

# gate (K values fixed; thresholds via --x0e/--x0f)
ACT_GATE_K_E = 15.0
ACT_GATE_K_F = 30.0

# STDP (verbatim)
TAU_PLUS, ALPHA, MU_PLUS, MU_MINUS = 20.0, 0.95, 0.4, 0.4
WMAX, WMAX_BS = 120.0, 30.0

# acceptance (k = 1.00 baseline)
BASELINE_CORR_F = -0.953


# ══════════════════════════════════════════════════════════════════════
# Metrics (pure numpy; also imported by summarize_small.py)
# ══════════════════════════════════════════════════════════════════════
def _corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = min(a.size, b.size)
    if n < 3:
        return float("nan")
    a, b = a[:n], b[:n]
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def resample_uniform(t_s, series):
    """Paced-gait sampling is non-uniform (tail chunks); resample like gate_diag."""
    t_s = np.asarray(t_s, float)
    d = np.diff(t_s)
    dt = float(np.median(d)) if d.size else 0.02
    if d.size and (d.max() - d.min()) / max(dt, 1e-12) > 0.05:
        tu = np.arange(t_s[0], t_s[-1] + 1e-12, dt)
        return tu, {k: np.interp(tu, t_s, np.asarray(v, float)) for k, v in series.items()}, dt
    return t_s, {k: np.asarray(v, float) for k, v in series.items()}, dt


def upcrossings(t, x):
    """Up-crossing times of x through its p10-p90 midpoint, with hysteresis
    (arm below 40 %, fire above 60 %), interpolated between samples."""
    x = np.asarray(x, float)
    if x.size < 8:
        return np.array([])
    lo, hi = np.percentile(x, 10), np.percentile(x, 90)
    if hi - lo < 1e-9:
        return np.array([])
    th_arm, th_fire, mid = lo + 0.4 * (hi - lo), lo + 0.6 * (hi - lo), 0.5 * (lo + hi)
    out, armed = [], False
    for i in range(1, x.size):
        if x[i] < th_arm:
            armed = True
        elif armed and x[i] >= th_fire:
            j = i
            while j > 0 and x[j - 1] >= mid:          # walk back to the mid crossing
                j -= 1
            if j > 0 and x[j] != x[j - 1]:
                f = (mid - x[j - 1]) / (x[j] - x[j - 1])
                out.append(t[j - 1] + f * (t[j] - t[j - 1]))
            else:
                out.append(t[j])
            armed = False
    return np.asarray(out, float)


def gate(r_hz, k, x0):
    z = np.clip(k * (np.asarray(r_hz, float) / RATE_REF_HZ - x0), -500, 500)
    return 1.0 / (1.0 + np.exp(-z))


def leg_metrics(t, tr, x0f, kf):
    fe, ff, rge, rgf = tr["force_e"], tr["force_f"], tr["rge"], tr["rgf"]
    p = lambda v, q: float(np.percentile(v, q))
    f10, f50, f90 = p(rgf, 10), p(rgf, 50), p(rgf, 90)
    span = f90 - f10
    up = upcrossings(t, fe)
    per = np.diff(up) if up.size > 2 else np.array([])
    n3 = max(1, rgf.size // 3)
    ampF_t1 = p(ff[:n3], 99) - p(ff[:n3], 1)
    ampF_t3 = p(ff[-n3:], 99) - p(ff[-n3:], 1)
    return dict(
        corr_F=_corr(fe, ff), corr_RG=_corr(rge, rgf),
        ampE=p(fe, 99) - p(fe, 1), ampF=p(ff, 99) - p(ff, 1),
        troughE=p(fe, 1), troughF=p(ff, 1), peakE=p(fe, 99), peakF=p(ff, 99),
        rge_p10=p(rge, 10), rge_p50=p(rge, 50), rge_p90=p(rge, 90),
        rgf_p10=f10, rgf_p50=f50, rgf_p90=f90, rgf_span=span,
        X0F_rec=0.5 * (f10 + f90) / RATE_REF_HZ,
        KF_rec=float(np.clip(600.0 / span, 3.0, 40.0)) if span > 1e-6 else float("nan"),
        df_p10=float(gate(f10, kf, x0f)), df_p90=float(gate(f90, kf, x0f)),
        df_swing=float(gate(f90, kf, x0f) - gate(f10, kf, x0f)),
        period_ms=float(per.mean() * 1000) if per.size else float("nan"),
        period_cv=float(per.std() / per.mean()) if per.size > 1 and per.mean() > 0 else float("nan"),
        n_cycles=int(per.size),
        rgf_shift_hz=p(rgf[-n3:], 50) - p(rgf[:n3], 50),
        ampF_ratio_t3_t1=float(ampF_t3 / ampF_t1) if ampF_t1 > 1e-6 else float("nan"),
        _upx=up,
    )


def lr_phase(up_l, up_r):
    """Phase of right-leg extensor onsets within the left-leg cycle, 0..1
    (0.5 = strict alternation). Circular mean and resultant length R."""
    if up_l.size < 3 or up_r.size < 3:
        return float("nan"), float("nan")
    T = float(np.mean(np.diff(up_l)))
    ph = []
    for tl in up_l[:-1]:
        nxt = up_r[up_r >= tl]
        if nxt.size:
            ph.append(((nxt[0] - tl) / T) % 1.0)
    if not ph:
        return float("nan"), float("nan")
    z = np.exp(2j * np.pi * np.asarray(ph))
    return float((np.angle(z.mean()) / (2 * np.pi)) % 1.0), float(abs(z.mean()))


def compute_metrics(times_ms, logs, x0f, kf, transient_s):
    """x0f: one threshold for both legs, or a dict {'L': .., 'R': ..}."""
    """logs[side][name] -> arrays. Returns {'leg_L':..., 'leg_R':..., 'lr_phase':...}."""
    t_s = np.asarray(times_ms, float) / 1000.0
    out = {}
    ups = {}
    for side in LEGS:
        names = ("force_e", "force_f", "rge", "rgf")
        tu, tr, _ = resample_uniform(t_s, {k: logs[side][k] for k in names})
        i0 = int(np.searchsorted(tu, tu[0] + transient_s))
        i0 = max(0, min(i0, tu.size - 8))
        x0 = x0f[side] if isinstance(x0f, dict) else x0f
        m = leg_metrics(tu[i0:], {k: v[i0:] for k, v in tr.items()}, x0, kf)
        ups[side] = m.pop("_upx")
        out[f"leg_{side}"] = m
    ph, R = lr_phase(ups["L"], ups["R"])
    out["lr_phase"] = ph
    out["lr_phase_R"] = R
    out["transient_s"] = float(transient_s)
    return out


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    return o


# ══════════════════════════════════════════════════════════════════════
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", default="40,40,20,20,20,40,40",
                    help="per-leg Izhikevich counts RGE,RGF,INE,INF,IAINT,MNE,MNF "
                         "(default = control m=0.40, 480 neurons both legs)")
    ap.add_argument("--relay-n", type=int, default=RELAY_N_DEFAULT,
                    help="size of each afferent/relay parrot population (not counted as neurons)")
    ap.add_argument("--conn-rule", default="bernoulli", choices=["bernoulli", "indegree"])
    ap.add_argument("--rgf-c", type=float, default=RGF_C, help="v1.3: RG-F Izhikevich c (reset, mV)")
    ap.add_argument("--rgf-d", type=float, default=RGF_D, help="v1.3: RG-F Izhikevich d (adaptation jump)")
    ap.add_argument("--i-e-rgf", type=float, default=I_E_RGF, help="v1.3: RG-F tonic current (pA)")
    ap.add_argument("--w-rg-rec-f", type=float, default=W_RG_REC_F, help="v1.3: RG-F recurrent weight")
    ap.add_argument("--indeg-min", type=int, default=1,
                    help="v1.2: minimum in-degree for Izhikevich-source projections (indegree rule)")
    ap.add_argument("--indeg-min-conserve", action="store_true",
                    help="v1.2: with --indeg-min, divide weights so mean input per neuron is unchanged")
    ap.add_argument("--preserve-input", action="store_true",
                    help="scale weights by N_ref/N_new so per-neuron input matches the control")
    ap.add_argument("--k-conn", type=float, default=0.20, help="connectivity scale (control 0.20)")
    ap.add_argument("--inh-comp", type=float, default=3.75,
                    help="gain on the E->InE->F pathway (RG-E->InE and InE->RG-F); 3.75 = k=0.20 recipe")
    ap.add_argument("--inh-comp-f", type=float, default=1.0,
                    help="gain on the F->InF->E pathway (RG-F->InF and InF->RG-E)")
    ap.add_argument("--w-ine", type=float, default=None,
                    help=f"InE->RG-F base weight before --inh-comp (default {W_INE2RGF})")
    ap.add_argument("--w-inf", type=float, default=None,
                    help=f"InF->RG-E base weight before --inh-comp-f (default {W_INF2RGE})")
    ap.add_argument("--stdp", default="off", choices=["off", "on"],
                    help="on = CUT->RG-E and BS->RG-E/F plastic at --stdp-lambda; "
                         "off = the same synapses and initial weights, lambda = 0")
    ap.add_argument("--stdp-lambda", type=float, default=1e-5)
    ap.add_argument("--sim-s", type=float, default=60.0)
    ap.add_argument("--transient-s", type=float, default=-1.0,
                    help="dropped before metrics; -1 = min(10, 20%% of --sim-s)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--x0f", type=float, default=1.42, help="flexor gate threshold (K_F fixed at 30)")
    ap.add_argument("--x0f-r", type=float, default=None,
                    help="right-leg flexor gate threshold (v1.1; default = --x0f)")
    ap.add_argument("--x0e", type=float, default=1.10, help="extensor gate threshold (K_E fixed at 15)")
    ap.add_argument("--out", default=".", help="output DIRECTORY")
    ap.add_argument("--tag", default="", help="file basename (default: derived from sizes/stdp/seed)")
    ap.add_argument("--cell", default="", help="cell label written to the JSON (grouping key)")
    ap.add_argument("--pass-label", default="", help="A or B (metadata)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--resolution-ms", type=float, default=0.2)
    ap.add_argument("--chunk-ms", type=float, default=20.0,
                    help="simulate chunk = trace sample interval (keep equal across compared runs)")
    ap.add_argument("--step-period-ms", type=float, default=520.0)
    ap.add_argument("--static-weight-cv", type=float, default=0.5)
    ap.add_argument("--delay-jitter-ms", type=float, default=0.2)
    ap.add_argument("--build-only", action="store_true",
                    help="build, print/save counts, exit without simulating")
    ap.add_argument("--nest-verbosity", default="M_ERROR")
    return ap


def main():
    args = build_parser().parse_args()
    t_wall0 = time.time()
    sizes = parse_sizes(args.sizes)
    k = float(args.k_conn)
    NR = int(args.relay_n)
    stdp_on = args.stdp == "on"
    lam = float(args.stdp_lambda) if stdp_on else 0.0
    SIM_MS = float(args.sim_s) * 1000.0
    X0F = {"L": float(args.x0f), "R": float(args.x0f if args.x0f_r is None else args.x0f_r)}
    TRANSIENT_S = (min(10.0, 0.2 * args.sim_s) if args.transient_s < 0 else float(args.transient_s))
    tag = args.tag or (f"s{'-'.join(str(sizes[x]) for x in SIZE_KEYS)}_{args.conn_rule}"
                       f"{'_pi' if args.preserve_input else ''}_stdp{args.stdp}_seed{args.seed}")
    os.makedirs(args.out, exist_ok=True)
    base = os.path.join(args.out, tag)

    # inhibition: direct overrides, then the pathway gains
    w_ine2rgf = float(args.w_ine) if args.w_ine is not None else W_INE2RGF
    w_inf2rge = float(args.w_inf) if args.w_inf is not None else W_INF2RGE
    g_e, g_f = float(args.inh_comp), float(args.inh_comp_f)
    W = dict(rg2ine=W_RG2INE * g_e, ine2rgf=w_ine2rgf * g_e,
             rg2inf=W_RG2INF * g_f, inf2rge=w_inf2rge * g_f)
    P = dict(rg_rec=P0_RG_REC * k, recip_f=P0_RG_RECIP_F * k, recip_e=P0_RG_RECIP_E * k,
             motor_recip=P0_MOTOR_RECIP * k, m2mus=P0_M2MUS * k, ia2in=P0_IA2IN * k,
             comm_f=P0_COMM_F * k, comm_e=P0_COMM_E * k)

    n_izh = 2 * izh_per_leg(sizes)
    print("=" * 72)
    print(f"{SCRIPT_VERSION}   tag={tag}   cell={args.cell or '-'}   pass={args.pass_label or '-'}")
    print(f"  sizes/leg : " + "  ".join(f"{kk}={sizes[kk]}" for kk in SIZE_KEYS)
          + f"   -> {izh_per_leg(sizes)}/leg, {n_izh} Izhikevich total (requested)")
    print(f"  relays    : {NR} per afferent/relay population (parrots, not counted)")
    print(f"  wiring    : k={k}  rule={args.conn_rule}  preserve_input={args.preserve_input}")
    print(f"  inhibition: W_RG2INE={W['rg2ine']:.2f} W_INE2RGF={W['ine2rgf']:.2f} "
          f"W_RG2INF={W['rg2inf']:.2f} W_INF2RGE={W['inf2rge']:.2f}")
    print(f"  STDP      : {args.stdp}  lambda={lam:g}")
    print(f"  gate      : X0_E={args.x0e} K_E={ACT_GATE_K_E}   X0_F L={X0F['L']} R={X0F['R']} "
          f"K_F={ACT_GATE_K_F}")
    print(f"  run       : {args.sim_s:g} s  seed={args.seed}  threads={args.threads}  "
          f"transient={TRANSIENT_S:g} s")
    print("=" * 72, flush=True)

    if nest is None:
        raise SystemExit("ERROR: NEST not importable (conda activate your NEST env)")
    if h5py is None and not args.build_only:
        raise SystemExit("ERROR: h5py not importable")
    try:
        nest.set_verbosity(args.nest_verbosity)
    except Exception:
        pass
    np.random.seed(args.seed)
    nest.ResetKernel()
    RES = float(args.resolution_ms)
    nest.SetKernelStatus({"resolution": RES, "local_num_threads": int(args.threads),
                          "print_time": False, "rng_seed": int(args.seed)})

    def q(x):
        return int(round(float(x) / RES)) * RES

    # delays (length_velocity, rat -- unchanged)
    def mk_delay(key, fallback):
        pr = DELAY_PRESETS["rat"].get(key)
        base_ms = (float(fallback) if pr is None else
                   pr["syn_delay_ms"] + pr["length_m"] / pr["velocity_mps"] * 1000.0)
        d = float(base_ms)
        if args.delay_jitter_ms > 0:
            d = d + nest.random.normal(mean=0.0, std=float(args.delay_jitter_ms))
        try:
            return nest.math.max(d, RES)
        except Exception:
            return max(float(base_ms), RES)
    DL = {kk: mk_delay(kk, fb) for kk, fb in [
        ("cut_to_rg", DELAY_MS), ("bs_to_rg", DELAY_MS), ("base_to_rg", DELAY_MS),
        ("rg_to_m", DELAY_MS), ("m_to_mus", DELAY_MS), ("ia_path", DELAY_MS),
        ("rg_rec", DELAY_MS), ("rg_recip", DELAY_RECIP_MS),
        ("motor_e2f", DELAY_MOTOR_RECIP_E2F_MS), ("motor_f2e", DELAY_MOTOR_RECIP_F2E_MS),
        ("commissural", DELAY_COMM_MS)]}

    # STDP initial weights: lognormal, CV 0.5, mean W0_IN (same for on and off)
    def winit(wmax):
        sigma = float(np.sqrt(np.log(1.0 + 0.5 ** 2)))
        mu = float(np.log(W0_IN) - 0.5 * sigma * sigma)
        p = nest.random.lognormal(mean=mu, std=sigma)
        return nest.math.min(nest.math.max(p, 0.0), float(wmax))
    W_INIT_CUT, W_INIT_BS = winit(WMAX), winit(WMAX_BS)

    #  build
    def relay(n, gen="poisson_generator"):
        g = nest.Create(gen, n)
        r = nest.Create("parrot_neuron", n)
        nest.Connect(g, r, conn_spec={"rule": "one_to_one"})
        return g, r

    leg = {}
    for side in LEGS:
        L = {}
        L["cut_pg"], L["cut_in"] = relay(NR)
        nest.SetStatus(L["cut_pg"], {"rate": CUT_RATE_OFF_HZ})
        L["bs_pg_e"], L["bs_in_e"] = relay(NR, "spike_generator")
        L["bs_pg_f"], L["bs_in_f"] = relay(NR, "spike_generator")
        period = 1000.0 / BS_REGULAR_HZ
        bt = np.arange(RES, SIM_MS + 1e-9, period)
        st = []
        for off in np.linspace(0.0, period, NR, endpoint=False):
            ta = bt + off
            ta = np.clip(ta + np.random.uniform(-BS_REGULAR_JITTER_MS, BS_REGULAR_JITTER_MS,
                                                size=ta.shape), RES, SIM_MS)
            ta = np.unique(np.sort(np.round(ta / RES) * RES))
            st.append(ta[ta >= RES].tolist())
        nest.SetStatus(L["bs_pg_e"], [{"spike_times": s} for s in st])
        nest.SetStatus(L["bs_pg_f"], [{"spike_times": s} for s in st])
        L["base_pg"], L["base_in"] = relay(NR)
        nest.SetStatus(L["base_pg"], {"rate": BASE_DRIVE_HZ})
        L["ia_pg_e"], L["ia_in_e"] = relay(NR)
        nest.SetStatus(L["ia_pg_e"], {"rate": IA_BASE_HZ_E})
        L["ia_pg_f"], L["ia_in_f"] = relay(NR)
        nest.SetStatus(L["ia_pg_f"], {"rate": IA_BASE_HZ_F})
        # stance Ia-E groups: as in the original, the POISSON GENERATORS project
        # to InE directly (independent trains per target; no parrot stage)
        L["ia_ext_pg_e"] = []
        for _ in range(3):
            g_ = nest.Create("poisson_generator", max(1, NR // 3))
            nest.SetStatus(g_, {"rate": 0.0})
            L["ia_ext_pg_e"].append(g_)
        L["ia_ext_pg_f"] = nest.Create("poisson_generator", NR)
        nest.SetStatus(L["ia_ext_pg_f"], {"rate": 0.0})

        L["rg_e"] = nest.Create("izhikevich", sizes["RGE"])
        L["rg_f"] = nest.Create("izhikevich", sizes["RGF"])
        L["m_e"] = nest.Create("izhikevich", sizes["MNE"])
        L["m_f"] = nest.Create("izhikevich", sizes["MNF"])
        L["ia_int_e"] = nest.Create("izhikevich", sizes["IAINT"])
        L["ia_int_f"] = nest.Create("izhikevich", sizes["IAINT"])
        L["in_e"] = nest.Create("izhikevich", sizes["INE"])
        L["in_f"] = nest.Create("izhikevich", sizes["INF"])
        for pop in ("rg_e", "rg_f", "m_e", "m_f"):
            nest.SetStatus(L[pop], izh_params)
        for pop in ("ia_int_e", "ia_int_f", "in_e", "in_f"):
            nest.SetStatus(L[pop], izh_inh_params)
        bias = LEFT_RIGHT_BIAS_IE if side == "L" else -LEFT_RIGHT_BIAS_IE
        nest.SetStatus(L["rg_e"], {"V_m": -65.0, "U_m": -13.0, "I_e": I_E_RGE})
        nest.SetStatus(L["rg_f"], {"a": RGF_A, "b": RGF_B, "c": float(args.rgf_c),
                                   "d": float(args.rgf_d), "V_m": -65.0, "U_m": RGF_B * -65.0,
                                   "I_e": float(args.i_e_rgf) + bias})
        nest.SetStatus(L["m_e"], {"V_m": -65.0, "U_m": -13.0, "I_e": I_E_MOTOR})
        nest.SetStatus(L["m_f"], {"V_m": -65.0, "U_m": -13.0, "I_e": I_E_MOTOR})
        L["mus_e"] = nest.Create("parrot_neuron", NR)
        L["mus_f"] = nest.Create("parrot_neuron", NR)
        for nm, pop in (("rge", "rg_e"), ("rgf", "rg_f"), ("ine", "in_e"), ("inf", "in_f"),
                        ("iainte", "ia_int_e"), ("iaintf", "ia_int_f"),
                        ("muse", "mus_e"), ("musf", "mus_f")):
            L["rec_" + nm] = nest.Create("spike_recorder")
            nest.Connect(L[pop], L["rec_" + nm])
        leg[side] = L

    stdp_base = {"tau_plus": TAU_PLUS, "lambda": lam, "alpha": ALPHA,
                 "mu_plus": MU_PLUS, "mu_minus": MU_MINUS, "Wmax": WMAX}
    for side in LEGS:
        nest.CopyModel("stdp_synapse", f"stdp_cut_rge_{side}", stdp_base)
        nest.CopyModel("stdp_synapse", f"stdp_bs_rge_{side}", {**stdp_base, "Wmax": WMAX_BS})
        nest.CopyModel("stdp_synapse", f"stdp_bs_rgf_{side}", {**stdp_base, "Wmax": WMAX_BS})

    # connect 
    # REF: izh population -> control size, used by --preserve-input.
    REF = {"rg_e": REF_SIZES["RGE"], "rg_f": REF_SIZES["RGF"], "in_e": REF_SIZES["INE"],
           "in_f": REF_SIZES["INF"], "ia_int_e": REF_SIZES["IAINT"],
           "ia_int_f": REF_SIZES["IAINT"], "m_e": REF_SIZES["MNE"], "m_f": REF_SIZES["MNF"]}
    projections = []      # bookkeeping for counts.json

    def C(name, side, src_key, tgt_key, p, w, d, *, src=None, tgt=None,
          model="static_synapse", rescale=True, device_src=False):
        """Connect with the chosen rule; rescale weight if --preserve-input and the
        source is a resized Izhikevich population. Returns expected in-degree."""
        Ls = leg[side]
        s = src if src is not None else Ls[src_key]
        t = tgt if tgt is not None else Ls[tgt_key]
        n_src = len(s)
        if p <= 0.0:
            return 0.0
        if args.conn_rule == "indegree" and not device_src:
            K_nat = min(n_src, max(1, int(round(p * n_src))))
            K = K_nat
            if args.indeg_min > 1 and rescale and src_key in REF:
                K = min(n_src, max(K_nat, int(args.indeg_min)))
            cs = {"rule": "fixed_indegree", "indegree": K, "allow_multapses": False}
            exp_in = float(K)
        else:
            cs = {"rule": "pairwise_bernoulli", "p": p}
            exp_in = p * n_src
        fac = 1.0
        n_ref = REF.get(src_key)
        if args.preserve_input and rescale and n_ref is not None and exp_in > 0:
            fac = (p * n_ref) / exp_in
        if args.conn_rule == "indegree" and not device_src and args.indeg_min_conserve \
                and K != K_nat:
            fac *= K_nat / K
        w_eff = w * fac
        nest.Connect(s, t, conn_spec=cs,
                     syn_spec={"synapse_model": model, "weight": w_eff, "delay": d})
        projections.append(dict(name=name, side=side, src=s, tgt=t, model=model,
                                rule=cs["rule"], p=p, exp_indegree=exp_in,
                                w_factor=fac,
                                w_mean=(float(w_eff) if isinstance(w_eff, (int, float)) else None)))
        return exp_in

    for side in LEGS:
        Ls = leg[side]
        C("CUT->RG-E", side, "cut_in", "rg_e", P_IN_STDP, W_INIT_CUT, DL["cut_to_rg"],
          model=f"stdp_cut_rge_{side}", rescale=False)
        C("BS->RG-E", side, "bs_in_e", "rg_e", P_IN_STDP, W_INIT_BS, DL["bs_to_rg"],
          model=f"stdp_bs_rge_{side}", rescale=False)
        C("BS->RG-F", side, "bs_in_f", "rg_f", P_IN_STDP, W_INIT_BS, DL["bs_to_rg"],
          model=f"stdp_bs_rgf_{side}", rescale=False)
        C("base->RG-E", side, "base_in", "rg_e", BASE_DRIVE_P, BASE_DRIVE_W, DL["base_to_rg"])
        C("base->RG-F", side, "base_in", "rg_f", BASE_DRIVE_P, BASE_DRIVE_W, DL["base_to_rg"])
        C("RG-E->M-E", side, "rg_e", "m_e", P_IN_STDP, W0_RM, DL["rg_to_m"])
        C("RG-F->M-F", side, "rg_f", "m_f", P_IN_STDP, W0_RM, DL["rg_to_m"])
        C("M-E->M-F", side, "m_e", "m_f", P["motor_recip"], W_MOTOR_RECIP, DL["motor_e2f"])
        C("M-F->M-E", side, "m_f", "m_e", P["motor_recip"], W_MOTOR_RECIP, DL["motor_f2e"])
        # parrot targets ignore weight; readout divides by the realised fan-in
        Ls["fanin_e"] = C("M-E->mus-E", side, "m_e", "mus_e", P["m2mus"], W_M2MUS,
                          DL["m_to_mus"], rescale=False)
        Ls["fanin_f"] = C("M-F->mus-F", side, "m_f", "mus_f", P["m2mus"], W_M2MUS,
                          DL["m_to_mus"], rescale=False)
        C("Ia-E->IaInt-E", side, "ia_in_e", "ia_int_e", IA2RG_P, W_IA_IN2INT, DL["ia_path"])
        C("IaInt-E->M-F", side, "ia_int_e", "m_f", IA2RG_P, W_IA_INT2ANT, DL["ia_path"])
        C("Ia-F->IaInt-F", side, "ia_in_f", "ia_int_f", IA2RG_P, W_IA_IN2INT, DL["ia_path"])
        C("IaInt-F->M-E", side, "ia_int_f", "m_e", IA2RG_P, W_IA_INT2ANT, DL["ia_path"])
        C("Ia-E->InE", side, "ia_in_e", "in_e", P["ia2in"], W_IA2IN, DL["ia_path"])
        C("Ia-F->InF", side, "ia_in_f", "in_f", P["ia2in"], W_IA2IN, DL["ia_path"])
        for gi, g_ in enumerate(Ls["ia_ext_pg_e"]):
            C(f"IaExt{gi}->InE", side, None, "in_e", P["ia2in"], W_IA2IN, DL["ia_path"],
              src=g_, device_src=True)
        C("flexAff->RG-F", side, "ia_ext_pg_f", "rg_f", P["ia2in"], W_FLEX_AFF2RGF,
          DL["ia_path"], device_src=True)
        C("flexAff->InF", side, "ia_ext_pg_f", "in_f", P["ia2in"], W_IA2IN,
          DL["ia_path"], device_src=True)
        C("RG-E->RG-E", side, "rg_e", "rg_e", P["rg_rec"], W_RG_REC_E, DL["rg_rec"])
        C("RG-F->RG-F", side, "rg_f", "rg_f", P["rg_rec"], float(args.w_rg_rec_f), DL["rg_rec"])
        C("RG-F->InF", side, "rg_f", "in_f", P["recip_f"], W["rg2inf"], DL["rg_recip"])
        C("InF->RG-E", side, "in_f", "rg_e", P["recip_f"], W["inf2rge"], DL["rg_recip"])
        C("RG-E->InE", side, "rg_e", "in_e", P["recip_e"], W["rg2ine"], DL["rg_recip"])
        C("InE->RG-F", side, "in_e", "rg_f", P["recip_e"], W["ine2rgf"], DL["rg_recip"])
        C("CUT->InE", side, "cut_in", "in_e", P_CUT2INE, W_CUT2INE, DL["cut_to_rg"])
    LL, RR = leg["L"], leg["R"]
    C("comm RG-F L->R", "L", "rg_f", None, P["comm_f"], W_COMM_F_INH, DL["commissural"], tgt=RR["rg_f"])
    C("comm RG-F R->L", "R", "rg_f", None, P["comm_f"], W_COMM_F_INH, DL["commissural"], tgt=LL["rg_f"])
    C("comm RG-E L->R", "L", "rg_e", None, P["comm_e"], W_COMM_E_INH, DL["commissural"], tgt=RR["rg_e"])
    C("comm RG-E R->L", "R", "rg_e", None, P["comm_e"], W_COMM_E_INH, DL["commissural"], tgt=LL["rg_e"])

    # lognormal heterogeneity on static weights (mean and sign preserved)
    scv = float(args.static_weight_cv)
    if scv > 0:
        sc = nest.GetConnections(synapse_model="static_synapse")
        if len(sc):
            w = np.asarray(sc.get("weight"), float)
            sig = float(np.sqrt(np.log(1 + scv * scv)))
            fac = np.random.default_rng(args.seed + 777).lognormal(-0.5 * sig * sig, sig, w.size)
            # one_to_one relay links (generator->parrot) have weight 1 and are
            # unaffected in effect: parrots repeat spikes regardless of weight.
            sc.set({"weight": (w * fac).tolist()})

    # ── ACTUAL counts from NEST ───────────────────────────────────────
    izh_ids = sorted(int(i) for side in LEGS for pop in
                     ("rg_e", "rg_f", "in_e", "in_f", "ia_int_e", "ia_int_f", "m_e", "m_f")
                     for i in leg[side][pop].tolist())
    izh_nc = nest.NodeCollection(izh_ids)
    try:
        n_izh_nest = len(nest.GetNodes({"model": "izhikevich"}))
    except Exception:
        n_izh_nest = -1
    per_pop = {side: {pop: len(leg[side][pop]) for pop in
                      ("rg_e", "rg_f", "in_e", "in_f", "ia_int_e", "ia_int_f", "m_e", "m_f")}
               for side in LEGS}
    syn_total = len(nest.GetConnections())
    syn_into_izh = len(nest.GetConnections(target=izh_nc))
    syn_izh_izh = len(nest.GetConnections(source=izh_nc, target=izh_nc))
    proj_counts = []
    for pr in projections:
        n = len(nest.GetConnections(source=pr["src"], target=pr["tgt"], synapse_model=pr["model"]))
        proj_counts.append(dict(name=pr["name"], side=pr["side"], n=n, rule=pr["rule"],
                                p=pr["p"], exp_indegree=pr["exp_indegree"],
                                w_factor=pr["w_factor"], n_tgt=len(pr["tgt"])))
    counts = dict(
        izh_total=len(izh_ids), izh_total_nest_getnodes=n_izh_nest,
        izh_per_leg={s: sum(per_pop[s].values()) for s in LEGS}, per_population=per_pop,
        synapses_total_all=syn_total,
        synapses_into_izh=syn_into_izh,
        synapses_izh_to_izh=syn_izh_izh,
        synapses_by_projection=proj_counts,
        relay_n=NR, sizes_requested=sizes)
    print(f"[counts] Izhikevich neurons built: {counts['izh_total']} "
          f"(GetNodes: {n_izh_nest})  per leg L={counts['izh_per_leg']['L']} "
          f"R={counts['izh_per_leg']['R']}")
    print(f"[counts] synapses: into Izhikevich={syn_into_izh}  Izh->Izh={syn_izh_izh}  "
          f"all (incl. relays/recorders)={syn_total}")
    for pc in proj_counts:
        if pc["side"] == "L" or pc["name"].startswith("comm"):
            print(f"   {pc['name']:<18}{pc['side']}  n={pc['n']:<6} rule={pc['rule']:<18} "
                  f"in/tgt={pc['n'] / max(1, pc['n_tgt']):6.2f}  w x{pc['w_factor']:.3f}")
    with open(base + ".counts.json", "w") as fh:
        json.dump(jsonable(counts), fh, indent=1)
    print(f"[counts] saved -> {base}.counts.json", flush=True)
    if args.build_only:
        return 0

    # simulate (per-chunk readout identical to cpg_optimization.py)
    CH = q(args.chunk_ms)
    keys = ("mus_e", "mus_f", "rge", "rgf", "ine", "inf", "iaint_e", "iaint_f",
            "act_e", "act_f", "force_e", "force_f", "len_e", "len_f", "ia_e", "ia_f")
    logs = {s: {kk: [] for kk in keys} for s in LEGS}
    state = {s: dict(act_e=0.0, act_f=0.0, force_e=0.0, force_f=0.0, len_e=L0, len_f=L0,
                     last={}) for s in LEGS}
    times = []
    stdp_models = [f"stdp_{x}_{s}" for s in LEGS for x in ("cut_rge", "bs_rge", "bs_rgf")]
    w_trace = {m: [] for m in stdp_models}
    w_times = []

    def cnt(side, nm):
        rec = leg[side]["rec_" + nm]
        cur = int(rec.get("n_events"))
        last = state[side]["last"].get(nm, 0)
        state[side]["last"][nm] = cur
        return cur - last

    def update_leg(side, dt_ms, cut_frac):
        Ls, S, Pl = leg[side], state[side], logs[side]
        dts = max(1e-9, dt_ms / 1000.0)
        r_me = cnt(side, "muse") / NR / dts / max(1.0, Ls["fanin_e"])
        r_mf = cnt(side, "musf") / NR / dts / max(1.0, Ls["fanin_f"])
        r_rge = cnt(side, "rge") / sizes["RGE"] / dts
        r_rgf = cnt(side, "rgf") / sizes["RGF"] / dts
        Pl["mus_e"].append(r_me); Pl["mus_f"].append(r_mf)
        Pl["rge"].append(r_rge); Pl["rgf"].append(r_rgf)
        Pl["ine"].append(cnt(side, "ine") / sizes["INE"] / dts)
        Pl["inf"].append(cnt(side, "inf") / sizes["INF"] / dts)
        Pl["iaint_e"].append(cnt(side, "iainte") / sizes["IAINT"] / dts)
        Pl["iaint_f"].append(cnt(side, "iaintf") / sizes["IAINT"] / dts)
        de = 1.0 / (1.0 + np.exp(-ACT_GATE_K_E * (r_rge / RATE_REF_HZ - args.x0e)))
        df = 1.0 / (1.0 + np.exp(-ACT_GATE_K_F * (r_rgf / RATE_REF_HZ - X0F[side])))
        ae = ACT_MAX * (1 - np.exp(-ACT_SAT_K * r_me)) * de
        af = ACT_MAX * (1 - np.exp(-ACT_SAT_K * r_mf)) * df
        ka = 1 - np.exp(-dts / (TAU_ACT_MS / 1000.0))
        S["act_e"] = min(ACT_MAX, max(0.0, S["act_e"] + ka * (ae - S["act_e"])))
        S["act_f"] = min(ACT_MAX, max(0.0, S["act_f"] + ka * (af - S["act_f"])))
        kF = 1 - np.exp(-dts / (TAU_FORCE_MS / 1000.0))
        tfe = FORCE_MAX * (1 - np.exp(-FORCE_SAT_K * S["act_e"]))
        tff = FORCE_MAX * (1 - np.exp(-FORCE_SAT_K * S["act_f"]))
        S["force_e"] = min(FORCE_MAX, max(0.0, S["force_e"] + kF * (tfe - S["force_e"])))
        S["force_f"] = min(FORCE_MAX, max(0.0, S["force_f"] + kF * (tff - S["force_f"])))
        kl = 1 - np.exp(-dts / (TAU_LENGTH_MS / 1000.0))
        S["len_e"] += kl * (L0 - S["len_e"]); S["len_f"] += kl * (L0 - S["len_f"])
        S["len_e"] -= SHORTEN_GAIN * S["force_e"] * dts
        S["len_f"] -= SHORTEN_GAIN * S["force_f"] * dts
        if cut_frac > 0:
            S["len_e"] += STRETCH_GAIN * cut_frac * dts
        S["len_e"] = min(L_MAX, max(L_MIN, S["len_e"]))
        S["len_f"] = min(L_MAX, max(L_MIN, S["len_f"]))
        ia_e = min(IA_RATE_MAX_HZ, max(0.0, IA_BASE_HZ_E + IA_K_FORCE_E * S["force_e"]
                                       + IA_K_STRETCH_E * max(0.0, S["len_e"] - L0)))
        ia_f = min(IA_RATE_MAX_HZ, max(0.0, IA_BASE_HZ_F + IA_K_FORCE_F * S["force_f"]
                                       + IA_K_STRETCH_F * max(0.0, S["len_f"] - L0)))
        Ls["ia_pg_e"].set(rate=ia_e); Ls["ia_pg_f"].set(rate=ia_f)
        for nm in ("act_e", "act_f", "force_e", "force_f", "len_e", "len_f"):
            Pl[nm].append(S[nm])
        Pl["ia_e"].append(ia_e); Pl["ia_f"].append(ia_f)

    clock = {"t": 0.0, "n": 0, "next_w": 0.0}

    def run_win(win, cut_frac):
        nc = int(win // CH)
        tail = q(win - nc * CH)
        for ci in range(nc + (1 if tail > 1e-9 else 0)):
            cur = q(CH if ci < nc else tail)
            if cur <= 0:
                continue
            nest.Simulate(cur)
            clock["t"] += cur; clock["n"] += 1
            for s in LEGS:
                update_leg(s, cur, cut_frac)
            times.append(clock["t"])
            if clock["t"] >= clock["next_w"]:          # weight means once per second
                w_times.append(clock["t"])
                for m in stdp_models:
                    w_trace[m].append(float(np.mean(nest.GetConnections(synapse_model=m).get("weight"))))
                clock["next_w"] += 1000.0

    SP = q(args.step_period_ms)
    HALF = q(SP / 2.0)
    SUB = q(HALF / 3.0)
    IEH = (60.0, 80.0, 100.0)
    IEF = 80.0
    NHC = max(2, int(np.ceil(SIM_MS / HALF)))
    t_sim0 = time.time()
    for hc in range(NHC):
        if clock["t"] >= SIM_MS - 1e-9:
            break
        stance, swing = ("L", "R") if hc % 2 == 0 else ("R", "L")
        leg[swing]["cut_pg"].set(rate=CUT_RATE_OFF_HZ)
        for g_ in leg[swing]["ia_ext_pg_e"]:
            g_.set(rate=0.0)
        leg[swing]["ia_ext_pg_f"].set(rate=IEF)
        leg[stance]["cut_pg"].set(rate=CUT_RATE_ON_HZ)
        leg[stance]["ia_ext_pg_f"].set(rate=0.0)
        for gi in range(3):
            leg[stance]["ia_ext_pg_e"][gi].set(rate=IEH[gi])
            if gi > 0:
                leg[stance]["ia_ext_pg_e"][gi - 1].set(rate=0.0)
            sr = min(SUB, SIM_MS - clock["t"])
            if sr <= 1e-9:
                break
            run_win(sr, 1.0)
        for g_ in leg[stance]["ia_ext_pg_e"]:
            g_.set(rate=0.0)
        if hc % 20 == 0:
            print(f"[sim] t={clock['t'] / 1000:.1f}/{args.sim_s:g} s  "
                  f"({time.time() - t_sim0:.0f} s wall)", flush=True)
    sim_wall = time.time() - t_sim0

    # metrics 
    met = compute_metrics(times, logs, X0F, ACT_GATE_K_F, TRANSIENT_S)
    wall = time.time() - t_wall0
    mL, mR = met["leg_L"], met["leg_R"]
    record = dict(
        script=SCRIPT_VERSION, cell=args.cell, pass_label=args.pass_label, tag=tag,
        created_utc=datetime.now(timezone.utc).isoformat(),
        nest_version=str(getattr(nest, "__version__", "?")), threads=int(args.threads),
        seed=int(args.seed), sim_s=float(args.sim_s), wall_s=wall, sim_wall_s=sim_wall,
        args=vars(args), sizes=sizes, izh_total=counts["izh_total"],
        synapses_into_izh=syn_into_izh, synapses_izh_to_izh=syn_izh_izh,
        counts=counts,
        gate=dict(X0_F=args.x0f, X0_F_R=X0F["R"], K_F=ACT_GATE_K_F, X0_E=args.x0e, K_E=ACT_GATE_K_E),
        # headline (left leg, as in FINDINGS)
        corr_F=mL["corr_F"], corr_RG=mL["corr_RG"],
        corr_F_R=mR["corr_F"], corr_RG_R=mR["corr_RG"],
        troughF=mL["troughF"], troughE=mL["troughE"], ampF=mL["ampF"], ampE=mL["ampE"],
        period_cv=mL["period_cv"], period_ms=mL["period_ms"],
        lr_phase=met["lr_phase"], lr_phase_R=met["lr_phase_R"],
        rgf_band=dict(p10=mL["rgf_p10"], p50=mL["rgf_p50"], p90=mL["rgf_p90"],
                      span=mL["rgf_span"], X0F_rec=mL["X0F_rec"], KF_rec=mL["KF_rec"],
                      df_swing=mL["df_swing"]),
        legs=met,
        stdp=dict(mode=args.stdp, lam=lam,
                  w_final_mean={m: (w_trace[m][-1] if w_trace[m] else None) for m in stdp_models}),
    )
    record["accept_run"] = dict(
        corr_F_ok=abs(mL["corr_F"] - BASELINE_CORR_F) <= 0.02,
        troughF_ok=mL["troughF"] <= 4.0,
        rg_gap_ok=abs(mL["corr_F"] - mL["corr_RG"]) <= 0.02)

    #  HDF5 (layout readable by gate_diag.py / plot_results.py) 
    with h5py.File(base + ".h5", "w") as h5:
        A = h5.attrs
        A["script"] = SCRIPT_VERSION
        A["cmdline"] = " ".join(sys.argv)
        A["nest_version"] = record["nest_version"]
        A["local_num_threads"] = int(args.threads); A["local_threads"] = int(args.threads)
        A["seed"] = int(args.seed)
        A["sim_ms"] = SIM_MS; A["dt_ms"] = CH; A["resolution_ms"] = RES
        A["transient_ms"] = TRANSIENT_S * 1000.0
        A["paced_gait"] = True; A["step_period_ms"] = SP
        A["k_connectivity"] = k; A["k_conn"] = k
        A["m_neurons"] = float(counts["izh_total"]) / 480.0 * 0.40   # equivalent m, info only
        A["izh_total"] = int(counts["izh_total"])
        A["sizes"] = args.sizes; A["conn_rule"] = args.conn_rule
        A["preserve_input"] = bool(args.preserve_input)
        A["INH_COMP"] = g_e; A["INH_COMP_F"] = g_f
        for kk, vv in W.items():
            A[f"W_{kk.upper()}"] = float(vv)
        A["stdp_mode"] = args.stdp; A["stdp_lambda"] = lam
        A["ACT_GATE_X0_F"] = args.x0f; A["ACT_GATE_K_F"] = ACT_GATE_K_F
        A["ACT_GATE_X0_F_R"] = X0F["R"]
        A["ACT_GATE_X0_E"] = args.x0e; A["ACT_GATE_K_E"] = ACT_GATE_K_E
        A["cell"] = args.cell; A["pass_label"] = args.pass_label
        A["counts_json"] = json.dumps(jsonable({kk: vv for kk, vv in counts.items()
                                                if kk != "synapses_by_projection"}))
        h5.create_dataset("times_ms", data=np.asarray(times, np.float32), compression="gzip")
        for side in LEGS:
            g = h5.create_group(f"leg_{side}")
            for kk, arr in logs[side].items():
                g.create_dataset(kk, data=np.asarray(arr, np.float32), compression="gzip")
        gw = h5.create_group("weights")
        gw.create_dataset("times_ms", data=np.asarray(w_times, np.float32))
        for m, v in w_trace.items():
            gw.create_dataset(m, data=np.asarray(v, np.float32))

    with open(base + ".json", "w") as fh:
        json.dump(jsonable(record), fh, indent=1)
    print(f"[done] {wall:.1f} s wall  corr_F(L)={mL['corr_F']:+.3f} corr_RG(L)={mL['corr_RG']:+.3f}  "
          f"corr_F(R)={mR['corr_F']:+.3f}  troughF={mL['troughF']:.2f}  ampF={mL['ampF']:.2f}  "
          f"L-R phase={met['lr_phase']:.3f}  RG-F p10/p90={mL['rgf_p10']:.1f}/{mL['rgf_p90']:.1f} Hz "
          f"-> X0F_rec L={mL['X0F_rec']:.3f} R={mR['X0F_rec']:.3f}")
    print(f"[out] {base}.h5  {base}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
