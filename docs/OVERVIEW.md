# Tiny two-leg spinal CPG — neuron and synapse reduction

Spiking CPG model of two legs (NEST 3.9, Izhikevich neurons). The question is
how far the network can be reduced — fewer synapses (`k`) and fewer neurons
(`m`) — while the flexor/extensor alternation still matches the full `k=1.00`
baseline, `corr_F = -0.953`.

**Result: 20 % of synapses can be removed. Neurons cannot.**
Recommended configuration: `k=0.20`, `m=0.40`, `inh_comp=3.75`, STDP `1e-5`.
480 neurons, 80 % of synapses, `corr_F = -0.94`.

See [`docs/RESULTS.md`](docs/RESULTS.md) for the numbers and
[`docs/METHOD.md`](docs/METHOD.md) for why every experiment is measured twice.

## Files

| file | role |
|---|---|
| `cpg_final_k025_m040_stdp.py` | the simulation |
| `cpg_adapter.py` | CLI shim in front of the model |
| `run_experiment.sh` | runs a grid, two passes |
| `gate_diag.py` | per-run metrics, prescribes each run's gate |
| `report.py` | per-condition table from a metrics CSV |
| `grids/*.tsv` | experiment definitions |

## Running an experiment

A grid is a TSV: condition label, then the flags that define it.

```
k200_m400	--k-conn 0.20 --m-neurons 0.40 --inh-comp 3.75 --stdp off
```

Everything else — timing, drive, delays, seeds — lives in `COMMON` inside
`run_experiment.sh` and is identical across every experiment.

```bash
conda activate nest39
bash run_experiment.sh grids/km.tsv check
bash run_experiment.sh grids/km.tsv passA
bash run_experiment.sh grids/km.tsv plan
bash run_experiment.sh grids/km.tsv passB
bash run_experiment.sh grids/km.tsv report
```

Output goes to `results_<grid>/passA` and `passB`; metrics to
`gate_<grid>_a.csv` and `gate_<grid>_b.csv`. Finished runs are skipped, so
stopping and resuming is free.

**Only pass B is evidence about the circuit.** Pass A exists to find where each
network's RG-F firing band sits; pass B re-runs each condition at the gate that
band warrants. See `docs/METHOD.md`.

## Settings

| variable | default | |
|---|---|---|
| `SEEDS` | 8 seeds | seed variance dominates below `m=0.40` |
| `SIM_MS` | 30000 | 30 s; use 60000–90000 with STDP |
| `ROOT` | `results_<grid>` | output directory |
| `FORCE` | 0 | 1 recomputes finished runs |
| `BASELINE` | -0.953 | the `k=1.00` reference |

## Environment

NEST 3.9.0, single-threaded. Thread count is part of the network definition —
NEST seeds its RNG per virtual process, so results are not comparable across
different thread counts or across x86/ARM.

```bash
conda create -n nest39 -c conda-forge nest-simulator=3.9 python=3.11 h5py matplotlib
```
