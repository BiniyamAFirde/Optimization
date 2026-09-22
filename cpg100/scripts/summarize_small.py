#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_small.py (v1.7.1) -- pass-B table, stability, symptoms, ladder, STDP settling and force traces
=====================================================================

  python3 summarize_small.py --root results_small                      (ladder, as v1.0)
  python3 summarize_small.py --root results_E0 --ref-prefix base_k100 \
          --targets base_k100=-0.953,ctrl_n480=-0.94                   (roadmap E0)

Reads <root>/passB/*.json (pass B ONLY -- pass A sits at the stock gate and is
not evidence), groups by cell, and reports mean +- SD over seeds.
"""
import argparse
import csv
import glob
import json
import os
import re
import statistics as st
import sys

import numpy as np

BASE_CORR = -0.953
BASE_ON = -0.973        # E0: k=1 baseline with STDP (1e-5), 60 s, 8 seeds, Mac NEST 3.9.0
A4_MODE = "one-sided"
C_E, C_F = "#d62728", "#1f77b4"          # plot1's colours


def load_runs(d):
    runs = []
    for p in sorted(glob.glob(os.path.join(d, "*.json"))):
        if p.endswith(".counts.json"):
            continue
        try:
            r = json.load(open(p))
        except Exception as e:
            print(f"[skip] {os.path.basename(p)}: {e!r}")
            continue
        r["_json"] = p
        r["_h5"] = p[:-5] + ".h5"
        runs.append(r)
    return runs


def num(x):
    return float("nan") if x is None else float(x)


def ms(v):
    v = [x for x in v if x == x]
    if not v:
        return float("nan"), float("nan")
    return st.mean(v), (st.stdev(v) if len(v) > 1 else 0.0)


def circ_mean(ph):
    ph = [p for p in ph if p == p]
    if not ph:
        return float("nan")
    z = np.mean(np.exp(2j * np.pi * np.asarray(ph)))
    return float((np.angle(z) / (2 * np.pi)) % 1.0)


def cell_order(c):
    m = re.search(r"n(\d+)", c)
    return (-(int(m.group(1)) if m else 0), "stdpon" in c, c)


def parse_targets(txt):
    out = {}
    for item in (txt or "").split(","):
        item = item.strip()
        if item:
            k, v = item.split("=")
            out[k.strip()] = float(v)
    return out


def third_stat(v, j, key):
    vals = []
    for x in v:
        th = (x.get("thirds") or {}).get("L")
        if th and len(th) == 3 and th[j] is not None:
            vals.append(num(th[j].get(key)))
    return ms(vals)[0]


def weight_drift(h5path):
    """max over plastic projections of |slope| over the last third, % per minute."""
    try:
        import h5py
        with h5py.File(h5path, "r") as h:
            if "weights" not in h:
                return float("nan")
            g = h["weights"]
            t = np.asarray(g["times_ms"][:], float) / 60000.0
            if t.size < 6:
                return float("nan")
            m = t >= t[0] + (t[-1] - t[0]) * 2 / 3
            out = []
            for k in g:
                if k == "times_ms":
                    continue
                w = np.asarray(g[k][:], float)
                if w.size != t.size or not np.isfinite(w[m]).all() or abs(w[m].mean()) < 1e-9:
                    continue
                slope = np.polyfit(t[m], w[m], 1)[0]
                out.append(abs(slope) / abs(w[m].mean()) * 100.0)
            return max(out) if out else float("nan")
    except Exception:
        return float("nan")


def summarize(runs, rg_criterion, ref_prefix, targets):
    by = {}
    for r in runs:
        by.setdefault(r.get("cell") or re.sub(r"_seed\d+$", "", r["tag"]), []).append(r)
    rows = []
    for cell in sorted(by, key=lambda c: (cell_order(c), c)):
        v = by[cell]
        g = lambda k: [num(x.get(k)) for x in v]
        cf, cf_sd = ms(g("corr_F"))
        crg, crg_sd = ms(g("corr_RG"))
        stab = [x.get("stability") for x in v]
        has_stab = all(isinstance(s, dict) and "corr_ok" in s for s in stab)
        n_stable = sum(1 for s in stab if isinstance(s, dict) and s.get("ok")) if has_stab else None
        n_cyc_ok = sum(1 for s in stab if isinstance(s, dict) and s.get("cycles_ok")) if has_stab else None
        d_corr_mean = ms([num(s.get("d_corrF_t3_t1")) for s in stab])[0] if has_stab else float("nan")
        row = dict(
            cell=cell, n=len(v), izh_total=v[0].get("izh_total"),
            syn_izh=int(round(ms(g("synapses_into_izh"))[0])),
            stdp=v[0].get("stdp", {}).get("mode", "?"),
            x0f=st.median(num(x["gate"]["X0_F"]) for x in v),
            corr_F=cf, corr_F_sd=cf_sd,
            d_base=cf - (BASE_ON if v[0].get("stdp", {}).get("mode") == "on" else BASE_CORR),
            corr_RG=crg, corr_RG_sd=crg_sd, gap=cf - crg,
            corr_F_R=ms(g("corr_F_R"))[0], corr_RG_R=ms(g("corr_RG_R"))[0],
            troughF=ms(g("troughF"))[0], troughF_sd=ms(g("troughF"))[1],
            troughE=ms(g("troughE"))[0], ampF=ms(g("ampF"))[0], ampE=ms(g("ampE"))[0],
            period_ms=ms(g("period_ms"))[0], period_cv=ms(g("period_cv"))[0],
            lr_phase=circ_mean(g("lr_phase")), lr_R=ms(g("lr_phase_R"))[0],
            rgf_p50=ms([num(x["rgf_band"]["p50"]) for x in v])[0],
            rgf_span=ms([num(x["rgf_band"]["span"]) for x in v])[0],
            cF_t1=third_stat(v, 0, "corr_F") if has_stab else float("nan"),
            cF_t2=third_stat(v, 1, "corr_F") if has_stab else float("nan"),
            cF_t3=third_stat(v, 2, "corr_F") if has_stab else float("nan"),
            ampF_t3_t1=(ms([num(s.get("ampF_t3_t1")) for s in stab])[0] if has_stab else float("nan")),
            rgf_p50_t1=third_stat(v, 0, "rgf_p50") if has_stab else float("nan"),
            rgf_p50_t3=third_stat(v, 2, "rgf_p50") if has_stab else float("nan"),
            n_stable=n_stable, n_cyc_ok=n_cyc_ok, d_corr_t3_t1=d_corr_mean,
            x0f_r=st.median(num(x["gate"].get("X0_F_R", x["gate"]["X0_F"])) for x in v),
            gap_R=ms(g("corr_F_R"))[0] - ms(g("corr_RG_R"))[0],
            wall_s=ms(g("wall_s"))[0],
            w_drift=(ms([weight_drift(x["_h5"]) for x in v])[0]
                     if v[0].get("stdp", {}).get("mode") == "on" else float("nan")),
            nest=",".join(sorted({str(x.get("nest_version")) for x in v})),
            threads=",".join(sorted({str(x.get("threads")) for x in v})),
            seeds=[x.get("seed") for x in v], _runs=v)
        rows.append(row)
    ref = {}
    for r in rows:
        if r["cell"].startswith(ref_prefix) and r["stdp"] not in ref:
            ref[r["stdp"]] = r
    for r in rows:
        c = ref.get(r["stdp"])
        r["ref_cell"] = c["cell"] if c else ""
        r["d_rg_ref"] = (r["corr_RG"] - c["corr_RG"]) if c else float("nan")
        a1 = abs(r["d_base"]) <= 0.02
        a2 = r["troughF"] <= 4.0
        a3 = r["corr_F_sd"] <= 0.03 and r["n"] > 1
        a4 = ((r["gap"] <= 0.02 if A4_MODE == "one-sided" else abs(r["gap"]) <= 0.02)
              if rg_criterion == "gap"
              else (abs(r["d_rg_ref"]) <= 0.02 if r["d_rg_ref"] == r["d_rg_ref"] else False))
        a5 = (None if r["n_stable"] is None else
              bool(r["d_corr_t3_t1"] <= 0.02 and r["ampF_t3_t1"] >= 0.90
                   and r["n_cyc_ok"] >= r["n"] - 1))
        r.update(A1=a1, A2=a2, A3=a3, A4=a4, A5=a5,
                 PASS=a1 and a2 and a3 and a4 and (a5 is not False))
        tgt = [(k, t) for k, t in targets.items() if r["cell"].startswith(k)]
        r["target"] = tgt[0][1] if tgt else float("nan")
        r["d_target"] = r["corr_F"] - r["target"] if tgt else float("nan")
        r["REPRO"] = (abs(r["d_target"]) <= 0.02) if tgt else None
        r["symptom"] = symptom(r)
    return rows


def symptom(r):
    """Roadmap symptom table -> next node. Empty for passing cells."""
    if r["PASS"]:
        return ""
    out = []
    if r["corr_RG"] > -0.80:
        out.append("no rhythm / weak alternation in RG (corr_RG %+.2f) -> E4 drive" % r["corr_RG"])
    elif r["gap"] > 0.02:
        out.append("readout off-centre: corr_F worse than corr_RG by %.3f -> re-gate" % r["gap"])
    if r["troughF"] > 4.0:
        out.append("shallow flexor trough %.2f -> E3 inhibition" % r["troughF"])
    if r["n"] > 1 and r["corr_F_sd"] > 0.03:
        out.append("seed SD %.3f -> too few inputs per neuron (fixed in-degree / raise p)" % r["corr_F_sd"])
    if r["A5"] is False:
        out.append("drifts over the run (dcorr %+.3f, ampF x%.2f) -> STDP settle / lower lambda"
                   % (r["d_corr_t3_t1"], r["ampF_t3_t1"]))
    if r.get("A8") is False:
        out.insert(0, "RHYTHM weaker than baseline at this gait (corr_RG %+.3f / R %+.3f, "
                      "baseline %+.3f worse by %.3f) -- force looks fine only because the gate "
                      "squares it -> E3/E5 at this gait" % (r["corr_RG"], r["corr_RG_R"],
                                                            r["corr_RG"] - r["d_rg_base"], r["d_rg_base"]))
    if out == [] and r.get("A6") is False:
        out.append("RIGHT leg off baseline (corr_F R %+.3f, corr_RG R %+.3f)" % (r["corr_F_R"], r["corr_RG_R"]))
    if out == [] and r.get("A7") is False:
        out.append("left-right phase off (%.3f, R=%.2f)" % (r["lr_phase"], r["lr_R"]))
    if not out and not r["A1"] and r["d_base"] > 0:
        out.append("alternation below baseline (corr_F %+.3f, corr_RG %+.3f) -> E3 inhibition"
                   % (r["corr_F"], r["corr_RG"]))
    if not out and not r["A1"] and r["d_base"] < 0:
        out.append("MORE alternating than baseline by %.3f (troughF %.2f): too much inhibition, "
                   "step it back" % (-r["d_base"], r["troughF"]))
    if not out and r["gap"] < -0.02:
        out.append("corr_F better than corr_RG by %.3f: steep gate squares the output (reference-type cell)"
                   % -r["gap"])
    return "; ".join(out) or "fails only by criterion flags shown"


def apply_gait_refs(rows):
    """E8: judge every cell against the k=1 STDP-off baseline at the same gait."""
    gait = lambda c: (int(m.group(1)) if (m := re.search(r"_g(\d+)(?:_|$)", c)) else None)
    base = {gait(r["cell"]): r for r in rows
            if r["cell"].startswith("base_k100") and r["stdp"] == "off" and gait(r["cell"])}
    shift = BASE_ON - BASE_CORR
    for r in rows:
        g = gait(r["cell"]); b = base.get(g)
        r["gait_ms"] = g
        if b is None:
            r["A6"] = r["A7"] = r["A8"] = None
            if not r["cell"].startswith("base_k100"):
                print(f"[warn] {r['cell']}: no _g<ms>_ in the name -> A6-A8 NOT checked, "
                      "A1 judged against the generic target", file=sys.stderr)
            continue
        extra = abs(shift) if r["stdp"] == "on" else 0.0
        inband = lambda x, ref: (ref - 0.02 - extra) <= x <= (ref + 0.02)
        r["target"], r["d_target"] = b["corr_F"], r["corr_F"] - b["corr_F"]
        r["d_base"] = r["corr_F"] - b["corr_F"]
        r["trough_lim"] = max(4.0, 1.25 * b["troughF"])
        r["A1"] = inband(r["corr_F"], b["corr_F"])
        r["A2"] = r["troughF"] <= r["trough_lim"]
        r["A6"] = bool(inband(r["corr_F_R"], b["corr_F_R"]) and r["gap_R"] <= 0.02)
        r["A7"] = bool(abs(r["lr_phase"] - 0.5) <= 0.05 and r["lr_R"] >= 0.90)
        r["d_rg_base"] = r["corr_RG"] - b["corr_RG"]
        r["d_rg_base_R"] = r["corr_RG_R"] - b["corr_RG_R"]
        r["A8"] = bool(r["d_rg_base"] <= 0.02 and r["d_rg_base_R"] <= 0.02)
        r["PASS"] = bool(r["A1"] and r["A2"] and r["A3"] and r["A4"]
                         and r["A5"] is not False and r["A6"] and r["A7"] and r["A8"])
        r["symptom"] = symptom(r) if not r["PASS"] else ""
    return base


def fmt_table(rows, rg_criterion):
    ok = lambda b: "n/a" if b is None else ("ok" if b else "--")
    hdr = ["cell", "n", "neur", "syn", "X0_F", "corr_F ± SD", "Δbase", "Δtarget", "corr_RG",
           "F−RG", "ΔRG ref", "corr_F R", "corr_RG R", "X0_F R", "troughF", "ampF", "periodCV", "L-R ph",
           "cF t1/t3", "ampF t3/t1", "stable", "w drift %/min", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "PASS"]
    lines = []
    for r in rows:
        stab = "n/a" if r["n_stable"] is None else f"{r['n_stable']}/{r['n']}"
        lines.append([r["cell"], str(r["n"]), str(r["izh_total"]), str(r["syn_izh"]),
                      f"{r['x0f']:.3f}", f"{r['corr_F']:+.3f} ± {r['corr_F_sd']:.3f}",
                      f"{r['d_base']:+.3f}",
                      ("" if r["REPRO"] is None else f"{r['d_target']:+.3f}"),
                      f"{r['corr_RG']:+.3f}", f"{r['gap']:+.3f}",
                      f"{r['d_rg_ref']:+.3f}", f"{r['corr_F_R']:+.3f}",
                      f"{r['corr_RG_R']:+.3f}", f"{r['x0f_r']:.3f}",
                      f"{r['troughF']:.2f}", f"{r['ampF']:.2f}", f"{r['period_cv']:.5f}",
                      f"{r['lr_phase']:.3f}",
                      ("n/a" if r["n_stable"] is None else f"{r['cF_t1']:+.3f}/{r['cF_t3']:+.3f}"),
                      ("n/a" if r["n_stable"] is None else f"{r['ampF_t3_t1']:.2f}"),
                      stab, ("" if r["w_drift"] != r["w_drift"] else f"{r['w_drift']:.2f}"),
                      ok(r["A1"]), ok(r["A2"]), ok(r["A3"]), ok(r["A4"]), ok(r["A5"]),
                      ok(r.get("A6")), ok(r.get("A7")), ok(r.get("A8")),
                      "PASS" if r["PASS"] else "fail"])
    w = [max(len(h), *(len(l[i]) for l in lines)) if lines else len(h) for i, h in enumerate(hdr)]
    txt = "  ".join(h.ljust(w[i]) for i, h in enumerate(hdr)) + "\n"
    txt += "  ".join("-" * w[i] for i in range(len(hdr))) + "\n"
    for l in lines:
        txt += "  ".join(l[i].ljust(w[i]) for i in range(len(hdr))) + "\n"
    md = "| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n"
    for l in lines:
        md += "| " + " | ".join(l) + " |\n"
    legend = (f"A1 |corr_F - baseline|<=0.02 (off -0.953, on {BASE_ON:+.3f})   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   "
              f"A4 {'corr_F-corr_RG<=+0.02' if (rg_criterion == 'gap' and A4_MODE == 'one-sided') else ('|corr_F-corr_RG|<=0.02' if rg_criterion == 'gap' else '|corr_RG - reference cell|<=0.02')}   "
              "A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)\n"
              "neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);\n"
              "all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);\n"
              "cF t1/t3 = corr_F in the first/last third of the post-transient window; "
              "L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).")
    return txt, md, legend


# ── traces ─────────────────────────────────────────────────────────────
def read_trace(path, leg="L"):
    import h5py
    with h5py.File(path, "r") as h:
        g = h[f"leg_{leg}"]
        fe = np.asarray(g["force_e"][:], float)
        ff = np.asarray(g["force_f"][:], float)
        if "times_ms" in h and len(h["times_ms"]) == fe.size:
            t = np.asarray(h["times_ms"][:], float) / 1000.0
        else:
            dt = float(h.attrs.get("dt_ms", 20.0)) / 1000.0
            t = np.arange(fe.size) * dt
    d = np.diff(t)
    dt = float(np.median(d))
    tu = np.arange(t[0], t[-1] + 1e-12, dt)
    return tu, np.interp(tu, t, fe), np.interp(tu, t, ff)


def window(t, fe, ff, t0, win):
    if t[-1] < t0 + win:
        t0 = max(t[0], t[-1] - win)
    m = (t >= t0) & (t < t0 + win)
    return t[m] - t0, fe[m], ff[m]


def representative(v):
    v = [x for x in v if os.path.exists(x["_h5"]) and x.get("corr_F") is not None]
    if not v:
        return None
    med = st.median(x["corr_F"] for x in v)
    return min(v, key=lambda x: abs(x["corr_F"] - med))


def pick_trace_cells(rows, wanted):
    if wanted:
        names = [w.strip() for w in wanted.split(",") if w.strip()]
        return [r for n in names for r in rows if r["cell"] == n]
    real = [r for r in rows if not r["cell"].startswith("chk")]
    out = [r for r in real if r["cell"].startswith(("ctrl", "base"))]
    out += [r for r in real if r["PASS"] and r not in out]
    seen = {r["izh_total"] for r in out if r["PASS"]}
    for r in real:
        if r not in out and not r["PASS"] and r["izh_total"] not in seen:
            out.append(r); seen.add(r["izh_total"])
    return out[:7]


def make_figure(rows, baseline, path, t0, win, wanted=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cells = pick_trace_cells(rows, wanted)
    ref = None
    if baseline and os.path.exists(baseline):
        t, fe, ff = read_trace(baseline)
        tb, feb, ffb = window(t, fe, ff, t0, win)
        m = t > min(10.0, 0.2 * t[-1])
        c_all = float(np.corrcoef(fe[m], ff[m])[0, 1]) if m.sum() > 10 else float("nan")
        ref = ("Baseline (k=1.00, 1200 neurons)", tb, feb, ffb, f"corr(Force-E,Force-F) = {c_all:+.3f}")
    nrow = len(cells) + (1 if ref else 0)
    fig, axes = plt.subplots(nrow, 1, figsize=(10.5, 2.25 * nrow), squeeze=False, sharey=True)
    axes = [a[0] for a in axes]
    ymax = 0.0
    i0 = 0
    if ref:
        ax = axes[0]
        name, tb, feb, ffb, lab = ref
        ax.plot(tb, feb, color=C_E, lw=1.5, label="Force-E")
        ax.plot(tb, ffb, color=C_F, lw=1.5, ls="--", label="Force-F")
        ax.set_title(f"{name}   {lab}", fontsize=9.5)
        ax.legend(loc="upper right", fontsize=7.5)
        ymax = max(ymax, feb.max(), ffb.max())
        i0 = 1
    for i, r in enumerate(cells):
        ax = axes[i0 + i]
        rep_ = representative(r["_runs"])
        if rep_ is None:
            ax.text(0.5, 0.5, "no .h5", ha="center", transform=ax.transAxes); continue
        t, fe, ff = read_trace(rep_["_h5"])
        tw, few, ffw = window(t, fe, ff, t0, win)
        ax.plot(tw, few, color=C_E, lw=1.5, label="Force-E")
        ax.plot(tw, ffw, color=C_F, lw=1.5, ls="--", label="Force-F")
        ax.axhline(4.0, color="#9aa0a6", lw=0.8, ls=":", zorder=0)
        ymax = max(ymax, few.max(), ffw.max())
        ax.set_title(f"{r['cell']}  ({r['izh_total']} neurons, {r['syn_izh']} syn, STDP {r['stdp']})   "
                     f"corr_F {r['corr_F']:+.3f} ± {r['corr_F_sd']:.3f}   corr_RG {r['corr_RG']:+.3f}   "
                     f"troughF {r['troughF']:.2f}   [{'PASS' if r['PASS'] else 'fail'}]   "
                     f"(seed {rep_.get('seed')}: {rep_['corr_F']:+.3f})", fontsize=8.3)
    for a in axes:
        a.set_ylim(0, max(18.0, ymax * 1.05)); a.set_ylabel("Force (a.u.)", fontsize=8)
        a.grid(True, color="#e6e6e6", lw=0.6)
    axes[-1].set_xlabel("time in window (s)")
    fig.suptitle("Muscle force alternation -- pass B (own gate); dotted = troughF limit 4.0; "
                 "one representative seed per cell (closest to the cell median corr_F)", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def make_stability_figure(rows, path):
    rows = [r for r in rows if r["n_stable"] is not None]
    if not rows:
        return False
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    cols = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ci, r in enumerate(rows):
        c = cols[ci % len(cols)]
        for k, x in enumerate(r["_runs"]):
            th = (x.get("thirds") or {}).get("L") or []
            if len(th) != 3 or any(t is None for t in th):
                continue
            xs = [1, 2, 3]
            lab = r["cell"] if k == 0 else None
            axes[0].plot(xs, [t["corr_F"] for t in th], "-o", color=c, alpha=0.55, ms=3, label=lab)
            axes[1].plot(xs, [t["ampF"] for t in th], "-o", color=c, alpha=0.55, ms=3)
            axes[2].plot(xs, [t["rgf_p50"] for t in th], "-o", color=c, alpha=0.55, ms=3)
    axes[0].axhline(BASE_CORR, color="#555", ls=":", lw=1)
    axes[0].axhspan(BASE_CORR - 0.02, BASE_CORR + 0.02, color="#999", alpha=0.15)
    for a, t in zip(axes, ("corr_F (leg L)", "ampF (leg L)", "RG-F median rate (Hz)")):
        a.set_xticks([1, 2, 3]); a.set_xticklabels(["1st third", "2nd", "3rd"])
        a.set_title(t, fontsize=10); a.grid(True, color="#e6e6e6", lw=0.6)
    axes[0].legend(fontsize=7.5, loc="lower left")
    fig.suptitle("Stability over the run, every seed (pass B). Flat lines = settled; "
                 "shaded = baseline ±0.02", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return True


def make_ladder_figure(rows, path):
    pts = {}
    for r in rows:
        m = re.search(r"(?:^|_)n(\d+)_(.+)$", r["cell"])
        if not m or r["cell"].startswith("chk"):
            continue
        var = ("ctrl " if r["cell"].startswith("ctrl") else "") + m.group(2)
        pts.setdefault(var, []).append((int(m.group(1)), r))
    if sum(len(v) for v in pts.values()) < 3:
        return False
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.0))
    for var in sorted(pts):
        v = sorted(pts[var])
        n = [x[0] for x in v]
        mk = "s" if var.startswith("ctrl") else "o"
        ax[0].errorbar(n, [x[1]["corr_F"] for x in v], yerr=[x[1]["corr_F_sd"] for x in v],
                       marker=mk, capsize=3, label=var)
        ax[1].plot(n, [x[1]["corr_RG"] for x in v], marker=mk, label=var)
        ax[2].plot(n, [x[1]["troughF"] for x in v], marker=mk, label=var)
    ax[0].axhspan(BASE_CORR - 0.02, BASE_CORR + 0.02, color="#999", alpha=0.18, lw=0)
    ax[0].axhline(BASE_CORR, color="#555", ls=":", lw=1)
    ax[1].axhspan(BASE_CORR - 0.02, BASE_CORR + 0.02, color="#999", alpha=0.18, lw=0)
    ax[2].axhline(4.0, color="#555", ls=":", lw=1)
    for a, t in zip(ax, ("corr_F (leg L), mean ± SD over seeds", "corr_RG (gate-independent rhythm)",
                         "troughF (flexor floor; limit 4.0)")):
        a.set_title(t, fontsize=10); a.set_xlabel("Izhikevich neurons (both legs)")
        a.invert_xaxis(); a.grid(True, color="#e6e6e6", lw=0.6)
    ax[0].legend(fontsize=7.5, loc="lower left")
    fig.suptitle("Size ladder, pass B. Shaded = baseline -0.953 ± 0.02", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    global A4_MODE, BASE_ON
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="results_small")
    ap.add_argument("--pass", dest="pas", default="B", choices=["A", "B"])
    ap.add_argument("--baseline", default="", help="k=1.00 baseline .h5 for a top trace panel")
    ap.add_argument("--rg-criterion", default="gap", choices=["gap", "control"])
    ap.add_argument("--a4", default="one-sided", choices=["one-sided", "symmetric"])
    ap.add_argument("--by-gait", action="store_true",
                    help="E8: judge cells ..._g<ms>_... against base_k100_g<ms>_stdpoff (adds A6, A7)")
    ap.add_argument("--base-on", type=float, default=BASE_ON,
                    help="baseline corr_F for STDP-on cells (E0: -0.973)")
    ap.add_argument("--ref-prefix", default="ctrl",
                    help="reference cell for ΔRG / --rg-criterion control (name prefix)")
    ap.add_argument("--targets", default="", help="reproduction targets, e.g. base_k100=-0.953,ctrl_n480=-0.94")
    ap.add_argument("--prefix", default="", help="output prefix (default <root>/small_summary)")
    ap.add_argument("--t0", type=float, default=20.0, help="start of the 5 s trace window (s)")
    ap.add_argument("--win", type=float, default=5.0)
    ap.add_argument("--trace-cells", default="", help="comma list of cells for the trace figure")
    args = ap.parse_args()

    A4_MODE, BASE_ON = args.a4, args.base_on
    d = os.path.join(args.root, f"pass{args.pas}")
    runs = load_runs(d)
    if not runs:
        sys.exit(f"no run JSONs in {d}")
    prefix = args.prefix or os.path.join(args.root, "small_summary" + ("" if args.pas == "B" else "_passA"))
    if args.pas == "A" and args.prefix:
        prefix += "_passA"
    targets = parse_targets(args.targets)
    rows = summarize(runs, args.rg_criterion, args.ref_prefix, targets)
    gait_base = apply_gait_refs(rows) if args.by_gait else {}

    nests = {v for r in rows for v in r["nest"].split(",")}
    thr = {v for r in rows for v in r["threads"].split(",")}
    if len(nests) > 1 or len(thr) > 1:
        print("!! MIXED PROVENANCE (NEST version or thread count differs) -- not comparable:",
              nests, thr)
    if nests - {"3.9.0"}:
        print(f"!! NEST {sorted(nests)} -- FINDINGS numbers were measured on 3.9.0; "
              "compare only within this run set.")

    txt, md, legend = fmt_table(rows, args.rg_criterion)
    print(f"\n=== pass {args.pas}: mean ± SD over seeds, baseline corr_F = {BASE_CORR} ===\n")
    print(txt)
    print(legend)
    if args.pas == "A":
        print("\n!! PASS A is at the stock gate -- read corr_RG and the band only.")

    verdict = []
    if targets:
        print("\n--- reproduction (pass B, leg L, tol 0.02) ---")
        for name, val in targets.items():
            cs = [r for r in rows if r["cell"].startswith(name)]
            if not cs:
                print(f"  {name:<12} target {val:+.3f}: NO CELLS FOUND")
                verdict.append(False)
                continue
            for r in cs:
                # v1.2: a reference fails only if the readout UNDER-reads the rhythm
                # (gap > +0.02); corr_F beating corr_RG is the steep gate squaring
                # a full-size burst, not a reproduction failure.
                good = (r["REPRO"] and r["A3"] and r["gap"] <= 0.02 and r["A5"] is not False)
                verdict.append(bool(good))
                why = []
                if not r["REPRO"]:
                    why.append(f"corr_F off target by {r['d_target']:+.3f}")
                if not r["A3"]:
                    why.append(f"only {r['n']} seed(s)" if r["n"] < 2
                               else f"seed SD {r['corr_F_sd']:.3f} > 0.03")
                if r["gap"] > 0.02:
                    why.append(f"corr_F-corr_RG gap {r['gap']:+.3f} (readout off-centre)")
                if r["A5"] is False:
                    why.append(f"drifts: mean dcorr_F {r['d_corr_t3_t1']:+.3f}, ampF x{r['ampF_t3_t1']:.2f}, "
                               f"cycles ok {r['n_cyc_ok']}/{r['n']}")
                print(f"  {r['cell']:<20} {r['corr_F']:+.3f} ± {r['corr_F_sd']:.3f} vs {val:+.3f}  "
                      f"corr_RG {r['corr_RG']:+.3f}  troughF {r['troughF']:.2f}  -> "
                      f"{'REPRODUCES' if good else 'does NOT reproduce: ' + '; '.join(why)}")
        print(f"\nE0 overall: {'REPRODUCES -> go to E1' if verdict and all(verdict) else 'does NOT reproduce -> fix pipeline'}")
    else:
        passing = [r for r in rows if r["PASS"]]
        if passing:
            best = min(passing, key=lambda r: r["izh_total"])
            print(f"\nSmallest passing cell: {best['cell']}  ({best['izh_total']} neurons, "
                  f"corr_F {best['corr_F']:+.3f} ± {best['corr_F_sd']:.3f})")
        else:
            print("\nNo cell passes all criteria.")

    drift = [r for r in rows if r["w_drift"] == r["w_drift"] and r["w_drift"] > 2.0]
    if drift:
        print("\n--- STDP settling (info) ---")
        for r in drift:
            print(f"  {r['cell']:<32} plastic weights still moving {r['w_drift']:.1f} %/min at the end "
                  f"(rhythm change t1->t3: {r['d_corr_t3_t1']:+.3f}) -> not settled; lower lambda or run longer")
    if gait_base:
        print("\n--- E8 verdict by gait (target = k=1 baseline at the same gait, +-0.02; STDP-on "
              f"cells may alternate up to {abs(BASE_ON - BASE_CORR):.3f} more, the E0 shift at 520 ms) ---")
        for g in sorted(gait_base):
            b = gait_base[g]
            print(f"  gait {g} ms: baseline corr_F {b['corr_F']:+.3f} (R {b['corr_F_R']:+.3f}), "
                  f"corr_RG {b['corr_RG']:+.3f} (R {b['corr_RG_R']:+.3f}), "
                  f"troughF {b['troughF']:.2f} -> trough limit {max(4.0, 1.25 * b['troughF']):.2f}")
            for r in rows:
                if r.get("gait_ms") == g and not r["cell"].startswith("base_k100"):
                    why = [k for k in ("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8") if r.get(k) is False]
                    print(f"     {r['cell']:<36} L {r['corr_F']:+.3f}±{r['corr_F_sd']:.3f}  "
                          f"R {r['corr_F_R']:+.3f}  RG {r['corr_RG']:+.3f}/{r['corr_RG_R']:+.3f}  "
                          f"trough {r['troughF']:.2f}  L-R {r['lr_phase']:.3f} "
                          f"(R={r['lr_R']:.2f})  -> {'PASS' if r['PASS'] else 'fail ' + ','.join(why)}")
        recipes = {}
        for r in rows:
            if r.get("gait_ms") and not r["cell"].startswith("base_k100"):
                key = re.sub(r"_g\d+", "", r["cell"])
                recipes.setdefault(key, []).append((r["gait_ms"], r["PASS"]))
        print("\n  recipe summary:")
        for k, v in recipes.items():
            ok_g = [g for g, p in sorted(v) if p]
            print(f"     {k:<34} passes at {ok_g if ok_g else 'no gait'} of {sorted(g for g, _ in v)}")
    fails = [r for r in rows if not r["PASS"] and not re.match(r"chk\d+_", r["cell"])
             and r["REPRO"] is None]
    if fails:
        print("\n--- symptoms (failing cells -> roadmap) ---")
        for r in fails:
            print(f"  {r['cell']:<24} {r['symptom']}")
    offR = [r for r in rows if r["gap_R"] > 0.02]
    if offR:
        print("\n--- right leg (info; pass rules use leg L until E8) ---")
        for r in offR:
            print(f"  {r['cell']:<24} corr_F R {r['corr_F_R']:+.3f} vs corr_RG R {r['corr_RG_R']:+.3f} "
                  f"(gap {r['gap_R']:+.3f}, X0_F R {r['x0f_r']:.3f}) -> right-leg readout off-centre")

    chk = [r for r in rows if re.match(r"chk\d+_", r["cell"])]
    if chk:
        print("\n--- readout-protocol checks (chkNN_<cell> vs <cell>) ---")
        byname = {r["cell"]: r for r in rows}
        for r in chk:
            twin = byname.get(re.sub(r"^chk\d+_", "", r["cell"]))
            if twin is None:
                print(f"  {r['cell']}: no matching cell"); continue
            print(f"  {r['cell']:<26} corr_F {r['corr_F']:+.3f}±{r['corr_F_sd']:.3f}  corr_RG {r['corr_RG']:+.3f}  "
                  f"troughF {r['troughF']:.2f}  ampF {r['ampF']:.2f}\n"
                  f"  {twin['cell']:<26} corr_F {twin['corr_F']:+.3f}±{twin['corr_F_sd']:.3f}  corr_RG {twin['corr_RG']:+.3f}  "
                  f"troughF {twin['troughF']:.2f}  ampF {twin['ampF']:.2f}\n"
                  f"  {'difference':<26} corr_F {r['corr_F'] - twin['corr_F']:+.3f}        corr_RG "
                  f"{r['corr_RG'] - twin['corr_RG']:+.3f}  troughF {r['troughF'] - twin['troughF']:+.2f}")

    with open(prefix + ".csv", "w", newline="") as fh:
        keys = [k for k in rows[0] if not k.startswith("_")]
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in keys})
    with open(prefix + ".md", "w") as fh:
        fh.write(f"# Summary -- pass {args.pas} ({args.root})\n\n{md}\n{legend}\n")
    print(f"\n[csv] {prefix}.csv\n[md]  {prefix}.md")
    try:
        make_figure(rows, args.baseline, prefix + "_force_traces.png", args.t0, args.win,
                    args.trace_cells)
        print(f"[png] {prefix}_force_traces.png")
    except Exception as e:
        print(f"[png] traces skipped: {e!r}")
    try:
        if make_ladder_figure(rows, prefix + "_ladder.png"):
            print(f"[png] {prefix}_ladder.png")
    except Exception as e:
        print(f"[png] ladder skipped: {e!r}")
    try:
        if make_stability_figure(rows, prefix + "_stability.png"):
            print(f"[png] {prefix}_stability.png")
    except Exception as e:
        print(f"[png] stability skipped: {e!r}")


if __name__ == "__main__":
    main()
