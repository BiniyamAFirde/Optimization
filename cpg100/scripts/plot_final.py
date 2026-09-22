#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_final.py -- plot1/plot2/plot3-style figures for one cell at 30 s and 60 s
=============================================================================

Rebuilds the project's three original figures (plot1 force, plot2 RG rates,
plot3 Ia rates) for a chosen cell, side by side with the k=1 baseline and the
480-neuron control, at TWO points of the run:

    left column   the 5 s window ending at 30 s  (25-30 s)
    right column  the 5 s window ending at 60 s  (55-60 s)

plus a full-run overview (0-60 s force, both legs, run thirds shaded) so drift
or a fading flexor would be visible. No new simulation: it reads the pass-B .h5
files the runners already wrote. For each cell it picks the seed whose corr_F is
closest to the cell's median (same rule as summarize_small.py), or --seed N.

  python3 plot_final.py                                  # defaults below
  python3 plot_final.py --cell n100_kmin3c_inex15_lam6_g520_stdpon --root results_E8 \
      --baseline results_E0 base_k100_stdpoff --control results_E6 ctrl_n480_stdpoff

Titles give corr(Force-E, Force-F) for the window shown AND for the whole run
after the 10 s transient (the number every table in the project uses).
Writes <out>_force.png, <out>_rg.png, <out>_ia.png, <out>_overview.png.
Needs numpy, h5py, matplotlib. No NEST.
"""
import argparse
import glob
import json
import os
import statistics as st
import sys

import numpy as np

C_E, C_F = "#d62728", "#1f77b4"      # plot1's colours (extensor red solid, flexor blue dashed)


def pick_run(root, cell, seed=None):
    js = [p for p in sorted(glob.glob(os.path.join(root, "passB", f"{cell}_seed*.json")))
          if not p.endswith(".counts.json")]
    runs = []
    for p in js:
        try:
            r = json.load(open(p))
        except Exception:
            continue
        if os.path.exists(p[:-5] + ".h5"):
            runs.append((r, p[:-5] + ".h5"))
    if not runs:
        return None
    if seed is not None:
        m = [x for x in runs if int(x[0].get("seed", -1)) == int(seed)]
        return m[0] if m else None
    med = st.median(r["corr_F"] for r, _ in runs)
    return min(runs, key=lambda x: abs(x[0]["corr_F"] - med))


def load(h5, leg="L"):
    import h5py
    with h5py.File(h5, "r") as h:
        t = np.asarray(h["times_ms"][:], float) / 1000.0
        g = h[f"leg_{leg}"]
        d = {k: np.asarray(g[k][:], float) for k in
             ("force_e", "force_f", "rge", "rgf", "ia_e", "ia_f") if k in g}
        step = float(h.attrs.get("step_period_ms", 520.0))
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1] + 1e-12, dt)
    return tu, {k: np.interp(tu, t, v) for k, v in d.items()}, step


def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1]) if a.size > 3 and a.std() > 0 and b.std() > 0 else float("nan")


def window(t, d, t_end, win):
    t_end = min(t_end, t[-1])
    m = (t > t_end - win) & (t <= t_end)
    return t[m] - (t_end - win), {k: v[m] for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell", default="n100_kmin3c_inex15_lam6_g520_stdpon")
    ap.add_argument("--root", default="results_E8")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--baseline", nargs=2, default=["results_E0", "base_k100_stdpoff"],
                    metavar=("ROOT", "CELL"), help="k=1 baseline (E0, 60 s, 520 ms)")
    ap.add_argument("--control", nargs=2, default=["results_E6", "ctrl_n480_stdpoff"],
                    metavar=("ROOT", "CELL"), help="480-neuron control; '- -' to omit")
    ap.add_argument("--ends", default="30,60", help="window end times (s)")
    ap.add_argument("--win", type=float, default=5.0)
    ap.add_argument("--leg", default="L", choices=["L", "R"])
    ap.add_argument("--out", default="", help="output prefix (default results_final/<cell>)")
    a = ap.parse_args()

    ends = [float(x) for x in a.ends.split(",")]
    rows = []
    for label, root, cell in (("Baseline k=1.00", *a.baseline), ("Control, 480 neurons", *a.control),
                              ("This work, 100 neurons", a.root, a.cell)):
        if root == "-":
            continue
        pr = pick_run(root, cell, a.seed if cell == a.cell else None)
        if pr is None:
            print(f"[skip] no pass-B run for {cell} in {root}")
            continue
        rec, h5 = pr
        t, d, step = load(h5, a.leg)
        m = t > min(10.0, 0.2 * t[-1])
        whole = corr(d["force_e"][m], d["force_f"][m])
        n_syn = rec.get("synapses_into_izh", "?")
        rows.append(dict(label=f"{label}  [{cell}, seed {rec.get('seed')}, {rec.get('izh_total')} neurons, "
                               f"{n_syn} syn, STDP {rec.get('stdp', {}).get('mode', '?')}, step {step:.0f} ms]",
                         short=label, t=t, d=d, whole=whole))
    if not rows or rows[-1]["short"] != "This work, 100 neurons":
        sys.exit(f"cell {a.cell} not found in {a.root}/passB -- check --root/--cell")

    out = a.out or os.path.join("results_final", a.cell)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    specs = [("force", ("force_e", "force_f"), ("Force-E", "Force-F"), "Force (a.u.)",
              "Muscle force alternation"),
             ("rg", ("rge", "rgf"), ("RG-E", "RG-F"), "RG rate (Hz)", "Rhythm-generator rates"),
             ("ia", ("ia_e", "ia_f"), ("Ia-E", "Ia-F"), "Ia rate (Hz)", "Ia afferent rates")]
    for tag, keys, labs, ylab, title in specs:
        if any(k not in rows[0]["d"] for k in keys):
            continue
        fig, axes = plt.subplots(len(rows), len(ends), figsize=(6.4 * len(ends), 2.6 * len(rows)),
                                 squeeze=False, sharey="row" if tag != "force" else True)
        for i, r in enumerate(rows):
            for j, te in enumerate(ends):
                ax = axes[i][j]
                tw, dw = window(r["t"], r["d"], te, a.win)
                ax.plot(tw, dw[keys[0]], color=C_E, lw=1.5, label=labs[0])
                ax.plot(tw, dw[keys[1]], color=C_F, lw=1.5, ls="--", label=labs[1])
                c = corr(dw[keys[0]], dw[keys[1]])
                ax.set_title(f"{r['short']}  |  {te - a.win:.0f}-{te:.0f} s   corr = {c:+.3f}"
                             + (f"   (whole run {r['whole']:+.3f})" if tag == "force" else ""),
                             fontsize=8.8)
                if j == 0:
                    ax.set_ylabel(ylab, fontsize=8.5)
                if i == len(rows) - 1:
                    ax.set_xlabel("time in window (s)")
                if tag == "force":
                    ax.axhline(4.0, color="#9aa0a6", lw=0.8, ls=":", zorder=0)
                ax.grid(True, color="#e6e6e6", lw=0.6)
                if i == 0 and j == len(ends) - 1:
                    ax.legend(loc="upper right", fontsize=7.5)
        if tag == "force":
            ymax = max(max(r["d"]["force_e"].max(), r["d"]["force_f"].max()) for r in rows)
            for a_ in fig.axes:
                a_.set_ylim(0, max(18.0, ymax * 1.05))
        fig.suptitle(f"{title} -- leg {a.leg}, 5 s windows at 30 s and 60 s (pass B, own gate)\n"
                     + "\n".join(r["label"] for r in rows), fontsize=8.5)
        fig.tight_layout()
        fig.savefig(f"{out}_{tag}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"[png] {out}_{tag}.png")

    # full-run overview of the chosen cell, both legs
    pr = pick_run(a.root, a.cell, a.seed)
    fig, axes = plt.subplots(2, 1, figsize=(15, 5.2), sharex=True)
    for ax, leg in zip(axes, ("L", "R")):
        t, d, step = load(pr[1], leg)
        ax.plot(t, d["force_e"], color=C_E, lw=0.7, label="Force-E")
        ax.plot(t, d["force_f"], color=C_F, lw=0.7, ls="--", label="Force-F")
        tr = min(10.0, 0.2 * t[-1])
        ax.axvspan(0, tr, color="#bbb", alpha=0.25, lw=0)
        for k in range(3):
            a0 = tr + k * (t[-1] - tr) / 3
            ax.axvline(a0, color="#888", lw=0.6, ls=":")
            m = (t >= a0) & (t < a0 + (t[-1] - tr) / 3)
            ax.text(a0 + 0.3, 17.2, f"third {k + 1}: corr {corr(d['force_e'][m], d['force_f'][m]):+.3f}",
                    fontsize=8)
        ax.set_ylim(0, 18.5); ax.set_ylabel(f"leg {leg}\nForce (a.u.)", fontsize=8.5)
        ax.grid(True, color="#e6e6e6", lw=0.6)
    axes[0].legend(loc="upper right", fontsize=7.5, ncol=2)
    axes[-1].set_xlabel("time (s)   (grey = 10 s transient excluded from all metrics)")
    fig.suptitle(f"{a.cell}, seed {pr[0].get('seed')}: the whole 60 s run, both legs", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{out}_overview.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[png] {out}_overview.png")


if __name__ == "__main__":
    main()
