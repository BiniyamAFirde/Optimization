#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import json
import os
import socket
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cpg_small as M  # noqa: E402  (the model; never edited)

E0_VERSION = "cpg_e0 1.1 (wraps cpg_small 1.0)"
EXPECTED_NEST = "3.9.0"

PRESETS = {
    "base_k100": ["--sizes", "100,100,50,50,50,100,100", "--relay-n", "100",
                  "--k-conn", "1.0", "--inh-comp", "1.0", "--inh-comp-f", "1.0",
                  "--conn-rule", "bernoulli", "--chunk-ms", "100"],
    "ctrl_n480": ["--sizes", "40,40,20,20,20,40,40", "--relay-n", "40",
                  "--k-conn", "0.20", "--inh-comp", "3.75", "--inh-comp-f", "1.0",
                  "--conn-rule", "bernoulli", "--chunk-ms", "100"],
}


def sha256(path):
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return None


def default_tag(a):
    s = M.parse_sizes(a.sizes)
    return (f"s{'-'.join(str(s[x]) for x in M.SIZE_KEYS)}_{a.conn_rule}"
            f"{'_pi' if a.preserve_input else ''}_stdp{a.stdp}_seed{a.seed}")


def thirds_metrics(h5_path, x0f, transient_s, step_period_ms):
    import h5py
    with h5py.File(h5_path, "r") as h:
        t_s = np.asarray(h["times_ms"][:], float) / 1000.0
        logs = {s: {k: np.asarray(h[f"leg_{s}"][k][:], float)
                    for k in ("force_e", "force_f", "rge", "rgf")} for s in M.LEGS}
    out = {}
    for side in M.LEGS:
        tu, tr, _ = M.resample_uniform(t_s, logs[side])
        i0 = int(np.searchsorted(tu, tu[0] + transient_s))
        n = tu.size - i0
        rows = []
        for j in range(3):
            a, b = i0 + j * n // 3, i0 + (j + 1) * n // 3
            if b - a < 8:
                rows.append(None)
                continue
            m = M.leg_metrics(tu[a:b], {k: v[a:b] for k, v in tr.items()}, x0f, M.ACT_GATE_K_F)
            m.pop("_upx", None)
            rows.append(dict(t0_s=float(tu[a]), t1_s=float(tu[b - 1]),
                             corr_F=m["corr_F"], corr_RG=m["corr_RG"], ampF=m["ampF"],
                             troughF=m["troughF"], troughE=m["troughE"],
                             rgf_p50=m["rgf_p50"], period_ms=m["period_ms"],
                             n_cycles=m["n_cycles"]))
        out[side] = rows
    L = out["L"]
    stab = dict(ok=False, reason="too short")
    if all(r is not None for r in L):
        dur = L[2]["t1_s"] - L[2]["t0_s"]
        exp_cyc = dur / (step_period_ms / 1000.0)
        d_corr = L[2]["corr_F"] - L[0]["corr_F"]
        amp_ratio = L[2]["ampF"] / L[0]["ampF"] if L[0]["ampF"] > 1e-9 else float("nan")
        cyc_ok = L[2]["n_cycles"] >= 0.9 * exp_cyc - 1
        stab = dict(d_corrF_t3_t1=d_corr, ampF_t3_t1=amp_ratio,
                    cycles_t3=L[2]["n_cycles"], cycles_expected_t3=exp_cyc,
                    corr_ok=bool(d_corr <= 0.02), amp_ok=bool(amp_ratio >= 0.90),
                    cycles_ok=bool(cyc_ok))
        stab["ok"] = bool(stab["corr_ok"] and stab["amp_ok"] and stab["cycles_ok"])
    return out, stab


def install_recorder_clearing():
    """Keep spike recorders from holding every event in RAM.

    cpg_small's readout only ever calls rec.get("n_events") on single spike
    recorders and differences successive values. This patch returns a running
    total kept here, and zeroes the recorder's in-memory buffer after each read
    (NEST: setting n_events = 0 clears the memory backend). The network, the spike
    counts and every trace are unchanged -- only the stored event lists go.
    Without it the k=1.00 baseline stores ~1.4e9 events (~20+ GB) over 60 s.
    """
    nest = M.nest
    NC = nest.NodeCollection
    if getattr(NC, "_e0_patched", False):
        return
    orig_get = NC.get
    cum = {}

    def get(self, *a, **k):
        if a == ("n_events",) and not k and len(self) == 1:
            cur = int(orig_get(self, "n_events"))
            gid = int(self.global_id)
            if cur:
                self.set(n_events=0)
                cum[gid] = cum.get(gid, 0) + cur
            return cum.get(gid, 0)
        return orig_get(self, *a, **k)

    NC.get = get            # one simulation per process, so gids never repeat
    NC._e0_patched = True


def main():
    # wrapper-only flags
    wp = __import__("argparse").ArgumentParser(add_help=False)
    wp.add_argument("--preset", choices=sorted(PRESETS), default=None)
    wp.add_argument("--strict-nest", action="store_true")
    wp.add_argument("--keep-recorder-events", action="store_true",
                    help="do NOT clear spike recorders (original memory behaviour)")
    wa, rest = wp.parse_known_args()
    if any(x in ("-h", "--help") for x in rest):
        print(__doc__)
        M.build_parser().print_help()
        return 0

    # preset first, explicit flags after -> explicit wins (argparse: last one wins)
    argv = (PRESETS[wa.preset] if wa.preset else []) + rest
    a = M.build_parser().parse_args(argv)

    nest_v = str(getattr(M.nest, "__version__", "?")) if M.nest is not None else "missing"
    if nest_v != EXPECTED_NEST:
        msg = (f"[E0] NEST {nest_v} != {EXPECTED_NEST}: results are NOT comparable with "
               f"the FINDINGS numbers (FINDINGS 7).")
        if wa.strict_nest:
            raise SystemExit(msg + " Aborting (--strict-nest).")
        print("!! " + msg, flush=True)

    if M.nest is not None and not wa.keep_recorder_events:
        install_recorder_clearing()
    sys.argv = [os.path.join(HERE, "cpg_small.py")] + argv
    rc = M.main()
    if rc not in (0, None) or a.build_only:
        return rc or 0

    tag = a.tag or default_tag(a)
    base = os.path.join(a.out, tag)
    transient = min(10.0, 0.2 * a.sim_s) if a.transient_s < 0 else float(a.transient_s)
    thirds, stab = thirds_metrics(base + ".h5", a.x0f, transient, a.step_period_ms)

    with open(base + ".json") as fh:
        rec = json.load(fh)
    rec["e0"] = dict(version=E0_VERSION, preset=wa.preset,
                     nest_version=nest_v, nest_expected=EXPECTED_NEST,
                     nest_ok=nest_v == EXPECTED_NEST, host=socket.gethostname(),
                     recorder_clearing=not wa.keep_recorder_events,
                     sha256_cpg_small=sha256(os.path.join(HERE, "cpg_small.py")),
                     sha256_cpg_e0=sha256(os.path.abspath(__file__)),
                     argv=argv)
    rec["thirds"] = thirds
    rec["stability"] = stab
    with open(base + ".json", "w") as fh:
        json.dump(M.jsonable(rec), fh, indent=1)
    L = thirds["L"]
    if all(r is not None for r in L):
        print("[E0] thirds (leg L)  corr_F " + " / ".join(f"{r['corr_F']:+.3f}" for r in L)
              + "   ampF " + " / ".join(f"{r['ampF']:.2f}" for r in L)
              + "   RG-F p50 " + " / ".join(f"{r['rgf_p50']:.1f}" for r in L)
              + f"   -> stable={'yes' if stab['ok'] else 'NO'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
