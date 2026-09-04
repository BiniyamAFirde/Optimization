#!/usr/bin/env python3

import sys, glob, os
import numpy as np
import h5py


def corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = min(len(a), len(b))
    if n < 3:
        return float("nan")
    a, b = a[:n], b[:n]
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def dom_freq(x, dt_s):
    x = np.asarray(x, float)
    if len(x) < 8:
        return float("nan")
    thr = x.mean()
    up = np.where((x[:-1] <= thr) & (x[1:] > thr))[0]
    if len(up) < 2:
        return float("nan")
    return float(1.0 / (np.mean(np.diff(up)) * dt_s))


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)

    files = []
    for pat in sys.argv[1:]:
        files.extend(glob.glob(pat))
    files = sorted(set(files))
    if not files:
        print("No files matched."); sys.exit(1)

    rows, empty = [], []
    for f in files:
        try:
            with h5py.File(f, "r") as h:
                if "leg_L" not in h:
                    empty.append(f)
                    continue
                g = h["leg_L"]
                fe = g["force_e"][:]; ff = g["force_f"][:]
                re_ = g["rge"][:];    rf = g["rgf"][:]
                dt_s = float(h.attrs.get("dt_ms", 100.0)) / 1000.0
                a = dict(h.attrs)
                rows.append(dict(
                    file=os.path.basename(f),
                    k=a.get("k_connectivity"), m=a.get("m_neurons"),
                    corrF=corr(fe, ff), corrRG=corr(re_, rf),
                    peakE=float(np.percentile(fe, 99)),
                    minE=float(np.percentile(fe, 1)),
                    minF=float(np.percentile(ff, 1)),
                    peakF=float(np.percentile(ff, 99)),
                    freq=dom_freq(fe, dt_s),
                    de_ptp=a.get("GATE_de_ptp"), df_ptp=a.get("GATE_df_ptp"),
                    rge_ptp=a.get("RG_RGE_ptp"), rgf_ptp=a.get("RG_RGF_ptp"),
                ))
        except Exception as e:
            print(f"[ERROR] {f}: {e!r}")

    if empty:
        print("\n!! These files contain NO datasets (the run's HDF5 write block aborted):")
        for f in empty:
            print("   ", f)
        print("   Re-run them; the metrics below exclude these.\n")

    if not rows:
        print("No usable files."); sys.exit(1)

    def fmt(v, w=8, p=3):
        return " " * w if v is None else f"{float(v):{w}.{p}f}"

    print(f"{'file':<52} {'corr(F)':>8} {'corr(RG)':>9} {'F_E lo':>7} {'F_E hi':>7} "
          f"{'F_F lo':>7} {'F_F hi':>7} {'freq':>6} {'de_ptp':>7} {'df_ptp':>7}")
    for r in rows:
        print(f"{r['file'][:52]:<52} {r['corrF']:8.3f} {r['corrRG']:9.3f} "
              f"{r['minE']:7.2f} {r['peakE']:7.2f} {r['minF']:7.2f} {r['peakF']:7.2f} "
              f"{r['freq']:6.3f}{fmt(r['de_ptp'],8)}{fmt(r['df_ptp'],8)}")

    print()
    for key, label in [("corrF", "corr(Force-E,Force-F)"),
                       ("corrRG", "corr(RG-E,RG-F)"),
                       ("peakE", "peak Force-E"),
                       ("minE", "trough Force-E  (baseline k=1.00 ~ 1.0)"),
                       ("minF", "trough Force-F  (baseline k=1.00 ~ 1.0)"),
                       ("freq", "frequency (Hz)")]:
        v = np.array([r[key] for r in rows], float)
        v = v[~np.isnan(v)]
        if v.size:
            print(f"{label:<42} min={v.min():8.4f}  max={v.max():8.4f}  mean={v.mean():8.4f}")


if __name__ == "__main__":
    main()
