#!/usr/bin/env python3

import argparse, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIGNALS = {
    "force": ("force_e", "force_f", "Force (a.u.)", "Muscle Force Alternation"),
    "rge":   ("rge", "rgf", "Rate (Hz)", "Rhythm-Generator (RG) Rates"),
    "act":   ("act_e", "act_f", "Activation", "Muscle Activation"),
    "ia":    ("ia_e", "ia_f", "Rate (Hz)", "Ia Sensory Generator Rates"),
}


def load(path, ke, kf):
    with h5py.File(path, "r") as h:
        if "leg_L" not in h:
            raise ValueError("file has no datasets (aborted write)")
        g = h["leg_L"]
        e = g[ke][:]; f = g[kf][:]
        dt_s = float(h.attrs.get("dt_ms", 100.0)) / 1000.0
        attrs = dict(h.attrs)
    return e, f, dt_s, attrs


def upsample(e, f, dt_s, factor):
    if factor <= 1:
        return e, f, dt_s
    n = len(e)
    x = np.arange(n)
    xi = np.linspace(0, n - 1, n * factor)
    try:
        from scipy.interpolate import CubicSpline
        return CubicSpline(x, e)(xi), CubicSpline(x, f)(xi), dt_s / factor
    except Exception:
        return np.interp(xi, x, e), np.interp(xi, x, f), dt_s / factor


def align(e, f):
    thr = (np.percentile(e, 1) + np.percentile(e, 99)) / 2.0
    up = np.where((e[:-1] <= thr) & (e[1:] > thr))[0]
    if len(up) == 0:
        return e, f
    return np.roll(e, -up[0]), np.roll(f, -up[0])


def label_for(path, attrs):
    base = os.path.basename(path)
    bits = []
    for a, nm in (("ACT_GATE_X0_E_used", "X0_E"), ("ACT_GATE_K_E_used", "K_E"),
                  ("ACT_GATE_X0_F_used", "X0_F"), ("ACT_GATE_K_F_used", "K_F"),
                  ("INH_COMP", "comp")):
        if a in attrs:
            bits.append(f"{nm}={float(attrs[a]):g}")
    return base[:46] + ("\n" + "  ".join(bits) if bits else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--baseline", default="")
    ap.add_argument("--signal", default="force", choices=list(SIGNALS))
    ap.add_argument("--window-s", type=float, default=5.0)
    ap.add_argument("--skip-s", type=float, default=10.0,
                    help="seconds to skip at the start (transient)")
    ap.add_argument("--smooth", type=int, default=1,
                    help="Upsample factor for DISPLAY ONLY via cubic interpolation. "
                         "Does not touch the simulation. Lowering the model's own "
                         "rate-update-ms to get smoother traces is NOT equivalent: the "
                         "rate feeds the activation gate and the Ia loop, so 20ms bins "
                         "made spike-count noise drive the dynamics and corr(RG) fell "
                         "from -0.91 to -0.55. Smooth the picture, not the model.")
    ap.add_argument("--out", default="force_runs.png")
    ap.add_argument("--overlay", action="store_true",
                    help="Overlay each run on the baseline in the SAME axes instead of "
                         "stacking separate panels. Baseline is drawn in grey behind the "
                         "run, so residual shape differences (peak height, floor depth, "
                         "burst width) are visible directly rather than by eye-matching "
                         "two rows. Traces are aligned on their first rising edge, since "
                         "runs have no common phase reference.")
    ap.add_argument("--metrics", action="store_true", default=True,
                    help="Annotate each panel with peak/floor for E and F vs baseline")
    args = ap.parse_args()

    ke, kf, ylab, title = SIGNALS[args.signal]

    files = []
    if args.baseline:
        b = sorted(glob.glob(args.baseline))
        if b:
            files.append((b[0], True))
    for pat in args.runs:
        for f in sorted(glob.glob(pat)):
            files.append((f, False))
    if not files:
        raise SystemExit("No files matched.")

    # overlay mode 
    if args.overlay:
        base = next((p for p, b in files if b), None)
        if base is None:
            raise SystemExit("--overlay needs --baseline")
        runs = [p for p, b in files if not b]
        be, bf, bdt, _ = load(base, ke, kf)
        i0 = int(args.skip_s / bdt); i1 = i0 + int(args.window_s / bdt)
        be, bf = be[i0:i1], bf[i0:i1]
        be, bf = align(be, bf)
        be, bf, bdt = upsample(be, bf, bdt, args.smooth)

        n = len(runs)
        fig, axes = plt.subplots(n, 1, figsize=(13, 3.0 * n), squeeze=False)
        fig.suptitle(f"{title} — overlaid on baseline (grey)", fontsize=14)
        for ax, path in zip(axes[:, 0], runs):
            try:
                e, f, dt_s, attrs = load(path, ke, kf)
            except Exception as ex:
                ax.text(0.5, 0.5, f"{os.path.basename(path)}\n{ex}", ha="center",
                        va="center", transform=ax.transAxes, color="crimson")
                continue
            j0 = int(args.skip_s / dt_s); j1 = j0 + int(args.window_s / dt_s)
            e, f = e[j0:j1], f[j0:j1]
            e, f = align(e, f)
            e, f, dt_s = upsample(e, f, dt_s, args.smooth)
            t = np.arange(len(e)) * dt_s
            tb = np.arange(len(be)) * bdt

            ax.plot(tb, be, color="0.65", lw=1.2, label="baseline E")
            ax.plot(tb, bf, color="0.65", lw=1.2, ls="--", label="baseline F")
            ax.plot(t, e, "r-", lw=1.6, label="E")
            ax.plot(t, f, "b--", lw=1.6, label="F")
            c = np.corrcoef(e, f)[0, 1] if e.std() > 1e-12 and f.std() > 1e-12 else np.nan
            note = (f"E {np.percentile(e,1):.1f}-{np.percentile(e,99):.1f}  "
                    f"F {np.percentile(f,1):.1f}-{np.percentile(f,99):.1f}   "
                    f"(base E {np.percentile(be,1):.1f}-{np.percentile(be,99):.1f}  "
                    f"F {np.percentile(bf,1):.1f}-{np.percentile(bf,99):.1f})")
            ax.set_title(f"{label_for(path, attrs)}   corr={c:.3f}\n{note}", fontsize=8)
            ax.set_ylabel(ylab, fontsize=8)
            ax.legend(fontsize=6, ncol=4, loc="upper right")
            ax.tick_params(labelsize=7)
        axes[-1, 0].set_xlabel("time in window (s)")
        fig.tight_layout(rect=[0, 0, 1, 0.98])
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        fig.savefig(args.out, dpi=130)
        print(f"saved -> {args.out}")
        return

    n = len(files)
    fig, axes = plt.subplots(n, 1, figsize=(12, 2.6 * n), squeeze=False)
    fig.suptitle(f"{title} — run comparison", fontsize=14)

    for ax, (path, is_base) in zip(axes[:, 0], files):
        try:
            e, f, dt_s, attrs = load(path, ke, kf)
        except Exception as ex:
            ax.text(0.5, 0.5, f"{os.path.basename(path)}\n{ex}",
                    ha="center", va="center", transform=ax.transAxes, color="crimson")
            ax.set_xticks([]); ax.set_yticks([])
            continue

        i0 = int(args.skip_s / dt_s)
        i1 = i0 + int(args.window_s / dt_s)
        e, f = e[i0:i1], f[i0:i1]
        t = np.arange(len(e)) * dt_s
        c = np.corrcoef(e, f)[0, 1] if e.std() > 1e-12 and f.std() > 1e-12 else float("nan")

        ax.plot(t, e, "r-", label=f"{ke}")
        ax.plot(t, f, "b--", label=f"{kf}")
        head = "BASELINE  " if is_base else ""
        ax.set_title(f"{head}{label_for(path, attrs)}   corr = {c:.3f}", fontsize=9)
        ax.set_ylabel(ylab, fontsize=8)
        ax.legend(fontsize=7, loc="upper right")
        ax.tick_params(labelsize=7)

    axes[-1, 0].set_xlabel("time in window (s)")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print(f"saved -> {args.out}")


if __name__ == "__main__":
    main()
