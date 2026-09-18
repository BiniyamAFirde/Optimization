#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpg_adapter.py -- run experiments 1-3 on the VALIDATED model
=============================================================

`run_experiments.sh` is written against `cpg_optimization.py`'s flag surface:
`--stdp off|cut_bs`, `--experiment`, `--run-label`, `--target-neurons`,
`--neuron-scope`, `--speed-label`, `--print-config`. Your working model,
`cpg_final_k025_m040_stdp.py`, has none of those — it uses `--stdp-lambda` and
`--m-neurons` — so pointing the driver at it fails with argparse exit 2 on every
run, as experiment 4 demonstrated.

This adapter presents `cpg_optimization.py`'s interface and forwards to the fork,
so:

    SCRIPT=$PWD/cpg_adapter.py bash run_experiments.sh local exp1

WHY NOT JUST USE cpg_optimization.py
------------------------------------
Because it has never been run. Every result in this study -- the k/m grid, the
gate prescription, the STDP ceiling -- came out of the fork. Swapping in an
unvalidated driver for the three headline experiments would make them
incomparable with everything already measured, which is the same mistake as
mixing thread counts or machines. The fork is the model of record.

WHAT IT TRANSLATES
------------------
  --stdp off                 -> --stdp-lambda 0
  --stdp cut_bs|cut_bs_ia|ia_only
                             -> --stdp-lambda <--stdp-lambda, default below>
                                (the fork's plastic set is already cut->RG-E and
                                bs->RG-E/F, which is exactly `cut_bs`; it has no
                                separate Ia-plastic mode, so cut_bs_ia and
                                ia_only are refused rather than silently
                                mistranslated)
  --target-neurons N --neuron-scope S
                             -> --m-neurons <m solved for N>, using
                                cpg_optimization.py's own population arithmetic
  --experiment/--run-label/--speed-label
                             -> dropped (metadata only; plot_results.py falls
                                back to the directory name and file stem, which
                                run_experiments.sh sets per experiment anyway)
  --print-config             -> resolved configuration printed, no simulation

Every other flag is passed through only if the fork's `--help` lists it, so an
unsupported extra is dropped with a note instead of killing the run.

STDP RATE -- READ THIS BEFORE RUNNING EXPERIMENT 2
--------------------------------------------------
`cpg_optimization.py` defaults `--stdp-lambda` to **1e-3**, and
`run_experiments.sh` does not override it, so "STDP ON" in experiment 2 means
lambda = 1e-3. Experiments 2d-2f measured that rate directly: at 1e-3 the flexor
amplitude falls to 0.79 and corr(F) to -0.647, and even at 1e-4 the network has
**not equilibrated after 90 s**. The usable ceiling was ~3e-6 at the stock gate
and ~1e-5 with the gate re-centred.

At 1e-3, experiment 2 will report that STDP destroys the rhythm -- but mostly
because the learning rate is about a hundred times too high, not because a small
network cannot carry plasticity. That confounds the question the experiment is
asking. This adapter therefore defaults to **1e-5** and says so in its banner.
To reproduce the brief's literal setting, pass CPG_STDP_LAMBDA=1e-3.
"""

import argparse
import os
import shlex
import subprocess
import sys

FORK = os.environ.get("CPG_MODEL", "cpg_final_k025_m040_stdp.py")
STDP_ON_LAMBDA = float(os.environ.get("CPG_STDP_LAMBDA", "1e-5"))


# ── cpg_optimization.py's population arithmetic, reproduced exactly ────
def population_sizes(m):
    def n(x):
        return max(1, int(x * m))
    rg_total = n(200)
    rg_e = max(1, rg_total // 2)
    return dict(N_CUT=n(100), N_BS=n(100), N_RG_TOTAL=rg_total, N_RG_E=rg_e,
                N_RG_F=max(1, rg_total - rg_e), N_MOTOR_E=n(100), N_MOTOR_F=n(100),
                N_MUS_E=n(100), N_MUS_F=n(100), N_IA_E=n(100), N_IA_F=n(100),
                N_IA_INT=n(50), N_INE=n(50), N_INF=n(50))


def neuron_counts(m, n_ia_groups=3, paced=True):
    s = population_sizes(m)
    izh_leg = (s["N_RG_E"] + s["N_RG_F"] + s["N_MOTOR_E"] + s["N_MOTOR_F"]
               + 2 * s["N_IA_INT"] + s["N_INE"] + s["N_INF"])
    relay_leg = (s["N_CUT"] + 2 * s["N_BS"] + s["N_BS"]
                 + s["N_IA_E"] + s["N_IA_F"] + s["N_MUS_E"] + s["N_MUS_F"])
    if paced:
        relay_leg += max(1, s["N_IA_E"] // max(1, n_ia_groups)) * max(1, n_ia_groups)
    out = dict(s)
    out.update(izh_per_leg=izh_leg, izh_total=2 * izh_leg,
               relay_total=2 * relay_leg, nodes_total=2 * (izh_leg + relay_leg),
               rg_total=2 * (s["N_RG_E"] + s["N_RG_F"]))
    return out


_SCOPE_KEY = {"izh_both_legs": "izh_total", "izh_per_leg": "izh_per_leg",
              "rg_only": "rg_total", "all_spiking": "nodes_total"}


def solve_m_for_target(target, scope, n_ia_groups=3, paced=True):
    key = _SCOPE_KEY[scope]
    best_m = best_err = best_n = None
    mi = 0.010
    while mi <= 1.0005:
        c = neuron_counts(round(mi, 4), n_ia_groups, paced)[key]
        err = abs(c - target)
        if best_err is None or err < best_err - 1e-12:
            best_m, best_err, best_n = round(mi, 4), err, c
        mi += 0.0005
    return best_m, int(best_n)


# ── CLI: cpg_optimization.py's surface ─────────────────────────────────
def build_parser():
    ap = argparse.ArgumentParser(
        description="Adapter: cpg_optimization.py's interface, the validated fork underneath")
    for f in ("--out", "--outdir", "--tag", "--run-name", "--experiment",
              "--run-label", "--speed-label", "--nest-verbosity", "--delay-model",
              "--species", "--save-weights", "--neuron-scope", "--stdp-winit-dist",
              "--sweep-pairs", "--sweep-dist"):
        ap.add_argument(f, default="")
    for f in ("--threads", "--seed", "--n-ia-groups", "--target-neurons",
              "--print-every", "--sweep-run-idx", "--max-weight-conns"):
        ap.add_argument(f, type=int, default=0)
    for f in ("--sim-ms", "--dt-ms", "--resolution-ms", "--simulate-chunk-ms",
              "--rate-update-ms", "--transient-ms", "--k-conn", "--m-neurons",
              "--inh-comp", "--inh-comp-f", "--w-ine2rgf", "--w-inf2rge",
              "--i-e-rgf", "--w-rg-rec-f", "--static-weight-cv", "--step-period-ms",
              "--stance-fraction", "--ia-ext-f-hz", "--ia-feedback-gain",
              "--cut-feedback-gain", "--act-gate-x0-e", "--act-gate-k-e",
              "--act-gate-x0-f", "--act-gate-k-f", "--delay-jitter-ms",
              "--delay-scale", "--bs-base-hz", "--bs-noise-std-hz",
              "--weight-sample-ms", "--stdp-alpha", "--stdp-tau-plus",
              "--stdp-mu-plus", "--stdp-mu-minus", "--wmax", "--wmax-bs",
              "--wmax-ia", "--p-ia2rg", "--stdp-winit-mean", "--stdp-winit-std",
              "--stdp-winit-min", "--stdp-winit-max"):
        ap.add_argument(f, type=float, default=None)
    ap.add_argument("--stdp-lambda", type=float, default=None)
    ap.add_argument("--ia-ext-hz", type=float, nargs="+", default=[])
    for f in ("--paced-gait", "--enforce-tonic-bs", "--long-run", "--print-config",
              "--weight-recorder", "--freeze-bs-rg", "--ablate-ia-loop",
              "--ablate-asym", "--ablate-comm", "--scale-taus-with-speed"):
        ap.add_argument(f, action="store_true")
    ap.add_argument("--stdp", default="off",
                    choices=["off", "cut_bs", "cut_bs_ia", "ia_only"])
    return ap


def main():
    args = build_parser().parse_args()

    if not os.path.isfile(FORK):
        sys.exit(f"[adapter] model not found: {FORK}\n"
                 f"          set CPG_MODEL=/path/to/cpg_final_k025_m040_stdp.py")
    help_txt = subprocess.run([sys.executable, FORK, "--help"],
                              capture_output=True, text=True).stdout

    # ── network size ───────────────────────────────────────────────────
    m = args.m_neurons if args.m_neurons is not None else 0.40
    solved = ""
    if args.target_neurons and args.target_neurons > 0:
        scope = args.neuron_scope or "izh_both_legs"
        if scope not in _SCOPE_KEY:
            sys.exit(f"[adapter] unknown --neuron-scope {scope}")
        m, got = solve_m_for_target(args.target_neurons, scope,
                                    args.n_ia_groups or 3, bool(args.paced_gait))
        solved = f"  (solved for ~{args.target_neurons} {scope} -> actual {got})"

    # ── plasticity ─────────────────────────────────────────────────────
    if args.stdp in ("cut_bs_ia", "ia_only"):
        sys.exit(f"[adapter] --stdp {args.stdp} has no equivalent in {FORK}: its "
                 f"plastic set is fixed at cut->RG-E and bs->RG-E/F. Use "
                 f"--stdp cut_bs, or run this condition on cpg_optimization.py.")
    if args.stdp == "off":
        lam = 0.0
    else:
        lam = args.stdp_lambda if args.stdp_lambda is not None else STDP_ON_LAMBDA

    counts = neuron_counts(m, args.n_ia_groups or 3, bool(args.paced_gait))
    if args.print_config or os.environ.get("CPG_VERBOSE_ADAPTER"):
        print("=" * 68)
        print(f"[adapter] model     : {FORK}")
        print(f"[adapter] network   : k={args.k_conn}  m={m}{solved}")
        print(f"[adapter]             RG-E/F {counts['N_RG_E']}/{counts['N_RG_F']} per leg, "
              f"{counts['izh_total']} spiking neurons both legs")
        print(f"[adapter] stdp      : {args.stdp} -> --stdp-lambda {lam:g}"
              + ("   (adapter default; see the module docstring)"
                 if args.stdp != "off" and args.stdp_lambda is None else ""))
        print(f"[adapter] duration  : {args.sim_ms} ms   step period {args.step_period_ms} ms")
        print(f"[adapter] out       : {args.out}")
        print("=" * 68)
    if args.print_config:
        return 0

    # ── build the fork's command line, dropping anything it lacks ──────
    cmd = [sys.executable, "-u", FORK]
    dropped = []

    def put(flag, *vals):
        if flag in help_txt:
            cmd.append(flag)
            cmd.extend(str(v) for v in vals)
        else:
            dropped.append(flag)

    put("--out", args.out)
    put("--m-neurons", m)
    put("--stdp-lambda", lam)
    if args.outdir:
        put("--outdir", args.outdir)
    if args.run_name:
        put("--run-name", args.run_name)
    if args.tag:
        put("--tag", args.tag)

    passthrough_num = [
        ("--k-conn", args.k_conn), ("--inh-comp", args.inh_comp),
        ("--inh-comp-f", args.inh_comp_f), ("--w-ine2rgf", args.w_ine2rgf),
        ("--w-inf2rge", args.w_inf2rge), ("--i-e-rgf", args.i_e_rgf),
        ("--w-rg-rec-f", args.w_rg_rec_f), ("--static-weight-cv", args.static_weight_cv),
        ("--sim-ms", args.sim_ms), ("--dt-ms", args.dt_ms),
        ("--resolution-ms", args.resolution_ms),
        ("--simulate-chunk-ms", args.simulate_chunk_ms),
        ("--rate-update-ms", args.rate_update_ms),
        ("--step-period-ms", args.step_period_ms),
        ("--stance-fraction", args.stance_fraction),
        ("--ia-ext-f-hz", args.ia_ext_f_hz),
        ("--act-gate-x0-e", args.act_gate_x0_e), ("--act-gate-k-e", args.act_gate_k_e),
        ("--act-gate-x0-f", args.act_gate_x0_f), ("--act-gate-k-f", args.act_gate_k_f),
        ("--delay-jitter-ms", args.delay_jitter_ms), ("--delay-scale", args.delay_scale),
        ("--bs-base-hz", args.bs_base_hz), ("--bs-noise-std-hz", args.bs_noise_std_hz),
        ("--weight-sample-ms", args.weight_sample_ms),
        ("--ia-feedback-gain", args.ia_feedback_gain),
        ("--cut-feedback-gain", args.cut_feedback_gain),
    ]
    for flag, v in passthrough_num:
        if v is not None:
            put(flag, v)
    for flag, v in (("--threads", args.threads), ("--seed", args.seed),
                    ("--n-ia-groups", args.n_ia_groups),
                    ("--max-weight-conns", args.max_weight_conns),
                    ("--print-every", args.print_every)):
        if v:
            put(flag, v)
    for flag, v in (("--nest-verbosity", args.nest_verbosity),
                    ("--delay-model", args.delay_model), ("--species", args.species),
                    ("--save-weights", args.save_weights)):
        if v:
            put(flag, v)
    if args.ia_ext_hz:
        put("--ia-ext-hz", *args.ia_ext_hz)
    for flag, on in (("--paced-gait", args.paced_gait),
                     ("--enforce-tonic-bs", args.enforce_tonic_bs),
                     ("--long-run", args.long_run),
                     ("--ablate-ia-loop", args.ablate_ia_loop),
                     ("--ablate-asym", args.ablate_asym),
                     ("--ablate-comm", args.ablate_comm)):
        if on:
            if flag in help_txt:
                cmd.append(flag)
            else:
                dropped.append(flag)

    if dropped:
        print(f"[adapter] not supported by {FORK}, dropped: {' '.join(sorted(set(dropped)))}")
    if os.environ.get("CPG_VERBOSE_ADAPTER"):
        print("[adapter] " + " ".join(shlex.quote(c) for c in cmd))

    return subprocess.call(cmd)


if __name__ == "__main__":
    sys.exit(main())
