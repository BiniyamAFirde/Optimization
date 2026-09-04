#!/usr/bin/env python3
"""
make_gate_table.py -- emit a per-seed gate table for the final runs.

"""
import sys, glob, re, os
import numpy as np
import h5py


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    files = []
    for pat in sys.argv[1:]:
        files.extend(glob.glob(pat))
    if not files:
        sys.stderr.write("No files matched.\n")
        sys.exit(1)

    rows = []
    for f in sorted(set(files)):
        try:
            with h5py.File(f, "r") as h:
                if "leg_L" not in h:
                    sys.stderr.write(f"skip (empty): {f}\n")
                    continue
                g = h["leg_L"]
                re_ = np.asarray(g["rge"][:], float)
                rf = np.asarray(g["rgf"][:], float)
                seed = int(dict(h.attrs).get("seed", 0))
        except Exception as e:
            sys.stderr.write(f"skip ({e!r}): {f}\n")
            continue

        if not seed:
            m = re.search(r"seed(\d+)", os.path.basename(f))
            seed = int(m.group(1)) if m else 0

        e1, e99 = np.percentile(re_, 1), np.percentile(re_, 99)
        f1, f99 = np.percentile(rf, 1), np.percentile(rf, 99)

     
        x0e = (e1 + 0.30 * (e99 - e1)) / 100.0
        x0f = (f1 + 0.50 * (f99 - f1)) / 100.0
        ke = float(np.clip(3.2 / max(((e99 - e1) / 2.0) / 100.0, 1e-6), 15, 20))
        kf = float(np.clip(3.2 / max(((f99 - f1) / 2.0) / 100.0, 1e-6), 28, 40))

        rows.append((seed, x0e, ke, x0f, kf, f1, f99))
        sys.stderr.write(
            f"seed {seed}: RG-E {e1:6.1f}-{e99:6.1f}  RG-F {f1:6.1f}-{f99:6.1f}"
            f"  -> x0e={x0e:.2f} ke={ke:.0f} x0f={x0f:.2f} kf={kf:.0f}\n")

    rows.sort()
    for seed, x0e, ke, x0f, kf, _, _ in rows:
        print(f"{seed}\t{x0e:.2f}\t{ke:.0f}\t{x0f:.2f}\t{kf:.0f}")

    if rows:
        f1s = np.array([r[5] for r in rows])
        f99s = np.array([r[6] for r in rows])
        sys.stderr.write(
            f"\nRG-F band varies across seeds: p1 spans {f1s.min():.1f}-{f1s.max():.1f} Hz, "
            f"p99 spans {f99s.min():.1f}-{f99s.max():.1f} Hz.\n"
            f"A single X0_F cannot sit at the midpoint of all of them -- that spread "
            f"is the size of the problem this table fixes.\n")


if __name__ == "__main__":
    main()
