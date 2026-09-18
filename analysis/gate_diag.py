#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
import h5py

RATE_REF_HZ = 100.0
SHIFT_BAD, SHIFT_WARN = 3.0, 1.8
RATIO_BAD, RATIO_WARN = 0.85, 0.93


def load(path, leg="L"):
    with h5py.File(path, "r") as h:
        if f"leg_{leg}" not in h:
            raise ValueError("no leg group")
        g = h[f"leg_{leg}"]
        names = ("force_e", "force_f", "rge", "rgf")
        tr = {k: np.asarray(g[k][:], float) for k in names if k in g}
        att = dict(h.attrs)
        t = None
        if "times_ms" in h:
            t = np.asarray(h["times_ms"][:], float) / 1000.0

    n = len(tr["force_e"])
    if t is None or len(t) != n:
        dt = float(att.get("dt_ms", 100.0)) / 1000.0
        t = np.arange(n) * dt

    def num(key, default=np.nan):
        v = att.get(key, default)
        return float(v) if v is not None else default

    meta = {
        "t": t,
        "sim_s": num("sim_ms", n * 100.0) / 1000.0,
        "X0_F": num("ACT_GATE_X0_F"),
        "K_F": num("ACT_GATE_K_F"),
        "lam": num("stdp_lambda", 0.0),
        "threads": int(num("local_num_threads", 0) or 0),
        "nest": str(att.get("nest_version", "?")),
    }
    return meta, tr, att


def gate(rate_hz, k, x0):
    z = k * (np.asarray(rate_hz, float) / RATE_REF_HZ - x0)
    return 1.0 / (1.0 + np.exp(-np.clip(z, -500.0, 500.0)))


def corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = min(len(a), len(b))
    if n < 3:
        return float("nan")
    a, b = a[:n], b[:n]
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def span(x):
    return float(np.percentile(x, 99) - np.percentile(x, 1))


def thirds(rgf, ff):
    n = len(rgf) // 3
    if n < 10:
        return [float("nan")] * 6
    out = []
    for i in range(3):
        s = slice(i * n, (i + 1) * n)
        out.append(float(np.median(rgf[s])))
    for i in range(3):
        s = slice(i * n, (i + 1) * n)
        out.append(span(ff[s]))
    return out


def settling(shift, amp1, amp3):
    ratio = amp3 / amp1 if amp1 > 1e-9 else 1.0
    if abs(shift) > SHIFT_BAD or ratio < RATIO_BAD:
        return "STILL DRIFTING"
    if abs(shift) > SHIFT_WARN or ratio < RATIO_WARN:
        return "marginal"
    return "settled"


def verdict(corr_rg, ampf_rel, span_rel):
    bad_rg = corr_rg > -0.80
    bad_read = ampf_rel < 0.75 and span_rel > 0.75
    if bad_rg and bad_read:
        return "MIXED  gate and circuit both degraded"
    if bad_rg:
        return "CIRCUIT  RG alternation lost"
    if bad_read:
        return "READOUT  gate mis-centred"
    return "OK  flexor preserved"


def analyse(meta, tr, transient_s=10.0):
    t = meta["t"]
    i0 = int(np.searchsorted(t, transient_s, "left"))
    i0 = min(i0, max(0, len(t) - 10))
    w = {k: v[i0:] for k, v in tr.items()}
    rge, rgf = w["rge"], w["rgf"]
    fe, ff = w["force_e"], w["force_f"]

    p10, p50, p90 = np.percentile(rgf, [10, 50, 90])
    rgf_span = float(p90 - p10)
    m1, m2, m3, a1, a2, a3 = thirds(rgf, ff)
    window_s = float(t[-1] - t[i0])
    shift = m3 - m1

    x0_rec = float(0.5 * (p10 + p90) / RATE_REF_HZ)
    k_rec = float(np.clip(600.0 / max(rgf_span, 1e-6), 3.0, 40.0))
    df = gate(rgf, meta["K_F"], meta["X0_F"])

    return {
        "lam": meta["lam"],
        "rgf_p10": float(p10), "rgf_p50": float(p50), "rgf_p90": float(p90),
        "rgf_span": rgf_span,
        "rgf_med_t1": m1, "rgf_med_t2": m2, "rgf_med_t3": m3,
        "ampF_t1": a1, "ampF_t2": a2, "ampF_t3": a3,
        "rgf_shift_hz": float(shift),
        "drift_window_s": window_s,
        "corr_RG": corr(rge, rgf),
        "corr_F": corr(fe, ff),
        "ampE": span(fe), "ampF": span(ff),
        "troughF": float(np.percentile(ff, 1)),
        "peakF": float(np.percentile(ff, 99)),
        "X0_F": meta["X0_F"], "K_F": meta["K_F"],
        "df_p10": float(np.percentile(df, 10)),
        "df_p50": float(np.percentile(df, 50)),
        "df_p90": float(np.percentile(df, 90)),
        "df_swing": float(np.percentile(df, 90) - np.percentile(df, 10)),
        "X0F_rec": x0_rec, "KF_rec": k_rec,
        "nest_version": meta["nest"], "threads": meta["threads"],
        "settling_state": settling(shift, a1, a3),
    }


FIELDS = [
    "file", "stem", "lam", "rgf_p10", "rgf_p50", "rgf_p90", "rgf_span",
    "rgf_med_t1", "rgf_med_t2", "rgf_med_t3", "ampF_t1", "ampF_t2", "ampF_t3",
    "rgf_shift_hz", "drift_window_s", "corr_RG", "corr_F", "ampE", "ampF",
    "troughF", "peakF", "X0_F", "K_F", "df_p10", "df_p50", "df_p90",
    "df_swing", "X0F_rec", "KF_rec", "nest_version", "threads",
    "span_rel", "ampF_rel", "verdict", "settling_state",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patterns", nargs="+")
    ap.add_argument("--leg", default="L")
    ap.add_argument("--transient-s", type=float, default=10.0)
    ap.add_argument("--csv", default="")
    args = ap.parse_args()

    files = []
    for pat in args.patterns:
        files.extend(glob.glob(pat))
    files = sorted(set(files))
    if not files:
        sys.exit("no files matched")

    rows = []
    for p in files:
        try:
            meta, tr, _ = load(p, args.leg)
            r = analyse(meta, tr, args.transient_s)
        except Exception as e:
            print(f"[skip] {os.path.basename(p)}: {e!r}")
            continue
        r["file"] = os.path.basename(p)
        r["stem"] = os.path.splitext(r["file"])[0]
        rows.append(r)

    if not rows:
        sys.exit("nothing readable")

    ref = rows[0]
    for r in rows:
        r["span_rel"] = r["rgf_span"] / max(ref["rgf_span"], 1e-9)
        r["ampF_rel"] = r["ampF"] / max(ref["ampF"], 1e-9)
        r["verdict"] = verdict(r["corr_RG"], r["ampF_rel"], r["span_rel"])

    hdr = f"{'run':<34}{'corr_RG':>9}{'corr_F':>9}{'ampF':>8}" \
          f"{'RG-F p50':>10}{'gate B':>9}  state"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['stem'][:34]:<34}{r['corr_RG']:>+9.3f}{r['corr_F']:>+9.3f}"
              f"{r['ampF']:>8.2f}{r['rgf_p50']:>10.1f}{r['X0F_rec']:>9.3f}"
              f"  {r['settling_state']}")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"\nwrote {args.csv}  ({len(rows)} runs)")


if __name__ == "__main__":
    main()
