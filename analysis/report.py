#!/usr/bin/env python3

import argparse
import csv
import os
import re
import statistics
import sys

KEYS = ("corr_F", "corr_RG", "ampF", "troughF", "rgf_span")


def group(path):
    if not os.path.exists(path):
        return {}
    by = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            label = re.sub(r"_seed\d+$", "", row["stem"])
            rec = {}
            for k in KEYS:
                try:
                    rec[k] = float(row[k])
                except (ValueError, KeyError):
                    rec[k] = float("nan")
            rec["state"] = row.get("settling_state", "")
            by.setdefault(label, []).append(rec)
    return by


def agg(runs, key):
    vals = [r[key] for r in runs if r[key] == r[key]]
    if not vals:
        return float("nan"), 0.0
    sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
    return statistics.mean(vals), sd


def states(runs):
    counts = {}
    for r in runs:
        counts[r["state"]] = counts.get(r["state"], 0) + 1
    return " ".join(f"{k}:{v}" for k, v in sorted(counts.items()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_file")
    ap.add_argument("--baseline", type=float, default=-0.953)
    args = ap.parse_args()

    by = group(args.csv_file)
    if not by:
        sys.exit(f"no rows in {args.csv_file}")

    print(f"\nbaseline k=1.00: corr_F = {args.baseline:+.3f}\n")
    hdr = f"{'condition':<26}{'n':>3}{'corr_F':>16}{'vs base':>9}" \
          f"{'corr_RG':>9}{'gap':>7}{'ampF':>7}{'troughF':>9}  state"
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for label in sorted(by):
        runs = by[label]
        cf, sd = agg(runs, "corr_F")
        crg, _ = agg(runs, "corr_RG")
        rows.append((label, len(runs), cf, sd, cf - args.baseline, crg,
                     cf - crg, agg(runs, "ampF")[0],
                     agg(runs, "troughF")[0], states(runs)))

    for lbl, n, cf, sd, d, crg, gap, amp, tro, st in rows:
        print(f"{lbl[:26]:<26}{n:>3}{cf:>+10.3f} +-{sd:>4.3f}{d:>+9.3f}"
              f"{crg:>+9.3f}{gap:>+7.3f}{amp:>7.2f}{tro:>9.2f}  {st}")

    best = min(rows, key=lambda r: r[2])
    print(f"\nbest: {best[0]}  corr_F {best[2]:+.3f} +-{best[3]:.3f}  "
          f"corr_RG {best[5]:+.3f}  troughF {best[8]:.2f}")

    near = [r for r in rows if r[4] < 0.03]
    if near:
        print("within 0.03 of baseline: " + ", ".join(r[0] for r in near))
    if abs(best[6]) > 0.05:
        print(f"gap {best[6]:+.3f}: readout is losing rhythm the circuit has")


if __name__ == "__main__":
    main()
