# cpg100 -- the two-leg spinal CPG reduced to 100 neurons

100 Izhikevich neurons (50 per leg) and 2,710 synapses (2.3 % of the k = 1 baseline's
1,200 neurons / 118,208 synapses) pass every acceptance rule at 350, 520 and 780 ms steps,
both legs, left-right phase, 60 s, 8 held-out seeds (STDP on). Report: `docs/CPG100_report.docx`;
details and limitations: `docs/FINDINGS_100neurons.md`; step-by-step log: `docs/ROADMAP_PROGRESS.md`.

| gait | baseline corr_F | 100 neurons corr_F |
|---|---|---|
| 350 ms | -0.918 +- 0.014 | -0.927 +- 0.019 |
| 520 ms | -0.955 +- 0.013 | -0.987 +- 0.005 |
| 780 ms | -0.994 +- 0.003 | -0.992 +- 0.003 |

## Reproduce (NEST 3.9.0, conda env `nest`)

    bash run_E8.sh --local all        # k=1 baselines at 3 gaits (slow: ~3 min per run)
    bash run_E8f.sh check
    bash run_E8f.sh --local all       # final recipe, 8 held-out seeds, ~30-40 min
    python3 plot_final.py --root results_E8f --cell n100_final_g520_stdpon \
        --baseline results_E8 base_k100_g520_stdpoff --control - -

Every `run_E*.sh` has `check`, `plan`, `smoke`, `--local` (laptop) and `submit` (Slurm arrays).
`cpg_small.py` is the model; `cpg_e*.py` are thin per-step drivers; the original
`cpg_2legs_fast.py` / `cpg_variantT_k025_m040.py` are unchanged.
Raw output (`results_*/`, `*.h5`) is not tracked; `results/` holds the summary tables and figures.
