# Tiny Two-Leg Spinal CPG — Network-Size Reduction

A spiking neural network model of a two-leg spinal central pattern generator (CPG),
reduced to a small fraction of its original neurons and synapses while preserving
flexor/extensor alternation. Built on [NEST](https://www.nest-simulator.org/).

## Result in one line

The circuit tolerates reduction to **k = 0.20 connectivity, m = 0.35 neuron count**
— 35 rhythm-generator neurons per half-centre, 17 interneurons, **~2.45% of the
baseline synapses** — while keeping extensor/flexor anti-correlation at
**corr(Force) = −0.957 ± 0.010** (8 seeds), against a full-network baseline of −0.953.

## Background

A spinal CPG produces rhythmic, alternating output for opposing muscle groups
(extensor vs flexor) via two reciprocally inhibiting rhythm-generator populations
(RG-E, RG-F). This project asks how far the network can be shrunk along two axes
before that alternation breaks:

- **k** — synaptic connectivity (scales all connection probabilities)
- **m** — neuron count (scales all population sizes)
- total synapses scale as **m²k**.

Two corrections make the reduced network work:

1. **Inhibitory weight compensation** (`--inh-comp`) — a *circuit-level* fix that
   scales up the surviving reciprocal-inhibition synapses to replace those lost to
   reduction (holds p·w roughly constant).
2. **Activation-gate recalibration** (`--act-gate-x0-e/-k-e/-x0-f/-k-f`) — a
   *readout-level* fix. The gate that converts RG firing rate into muscle force was
   calibrated for the intact network's firing range and mis-centres after reduction.
   Recalibrating it — per network for best results — recovers alternation.

Full write-up: `docs/CPG_reduction_report.pdf`.

## Repository layout

| Folder | Contents |
|---|---|
| `models/` | The CPG model scripts (see below) |
| `analysis/` | Metrics, gate calibration, and plotting |
| `optimize/` | CMA-ES automatic parameter search |
| `slurm/` | HPC batch scripts |
| `figures/` | Result figures |
| `docs/` | Report (PDF + Word) |

### Model scripts

| Script | Point | corr(Force) |
|---|---|---|
| `cpg_2legs_fast.py` | full network (baseline, STDP available) | −0.953 |
| `cpg_2legs_k020_m035.py` | **recommended** reduced point | −0.957 ± .010 |
| `cpg_2legs_k025_m040.py` | conservative reduced point | −0.969 ± .005 |
| `cpg_2legs_opt.py` | STDP-off model used by the optimizer | — |
| `cpg_2legs_reduced_stdp.py` | reduction + STDP (for future learning experiments) | — |

The two reduced scripts have the recipe baked into their defaults, so they
reproduce with no flags.

## Requirements

- Python 3.9+
- NEST 3.x (`pip` cannot install this; see the [NEST install guide](https://nest-simulator.readthedocs.io/en/stable/installation/index.html))
- `pip install -r requirements.txt` for the rest

## Quick start

**Run one reduced-network simulation** (single-threaded — see note):

```bash
python3 models/cpg_2legs_k020_m035.py \
    --out cpg_run.h5 --outdir results/ \
    --seed 12345 --threads 1 \
    --sweep-pairs "3.5:0.30" --sweep-run-idx 0 --sweep-dist lognormal_cv \
    --sim-ms 30000 --dt-ms 10 --paced-gait --step-period-ms 520 \
    --save-weights none --long-run
```

**Score it:**

```bash
python3 analysis/extract_metrics4.py "results/*.h5"
```

**Plot it against baseline** (`--smooth` is display-only interpolation):

```bash
python3 analysis/plot_force_runs.py --overlay --smooth 8 \
    --baseline "results_baseline/*.h5" \
    --runs "results/*.h5" --out figures/overlay.png
```

### Per-network gates (the best result)

The −0.957 figure uses gates calibrated *per network* from measured firing rates:

```bash
# 1. run 8 seeds with the fixed default gate
sbatch slurm/run_k020_m035.sh
# 2. derive each network's own gate from its RG firing band
python3 analysis/make_gate_table.py "results_stage8/*.h5" > gates.tsv
# 3. re-run with per-network gates (script reads gates.tsv)
sbatch slurm/run_k020_m035.sh    # per-network variant
```

## Important: run single-threaded

NEST is **not reproducible at >1 thread** for these small populations — identical
settings and seed can differ by up to 0.08 in corr(Force). Always pass
`--threads 1` (and `--nest-seed` where the flag exists). This is a known,
documented behaviour, not a bug in this code.

## Reproducing the figures

`figures/` were produced from the 8-seed per-network run (`results_stage9/`):

- `s9_force_smooth.png` — force alternation, all 8 seeds vs baseline
- `s9_rge_best.png` — rhythm-generator rates (shows RG-F goes tonic — the key mechanism finding)
- `s9_ia_best.png` — Ia sensory rates

## Automatic parameter search (optional)

NEST is not differentiable, so this uses **CMA-ES** (gradient-free), not gradient
descent. The objective is *constrained* — it rewards corr(Force) only while the
rhythm and force floors stay healthy, so it cannot fake alternation from a dead
rhythm.

```bash
sbatch optimize/run_optimize.sh    # searches at k=0.25/m=0.40
```

## Key findings

- **The flexor fails first.** Every failure mode traces to RG-F, never RG-E.
- **The mechanism changes.** At reduced k, RG-F stops bursting and goes tonic;
  alternation is recovered by compensated inhibition + rescaled readout, not by a
  restored half-centre.
- **Reduction floor:** k < 0.15 or m < 0.30, where RG-F's firing range collapses
  below ~15 Hz and no gate setting recovers alternation.
- **corr(Force) alone can mislead** — at k=0.10 the gate faked −0.95 from a broken
  rhythm (corr(RG) −0.78). Always read corr(RG) alongside it.

## Citing / contact

Internship project — see `docs/` for the full report and authors.
