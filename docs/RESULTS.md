# Results

NEST 3.9.0, 1 thread, paced gait 520 ms, 8 seeds, 30 s per run unless stated.
Every condition measured at its own gate (`docs/METHOD.md`).
Baseline: `k=1.00`, `corr_F = -0.953`.

## Recommended configuration

**`k = 0.20`, `m = 0.40`, `inh_comp = 3.75`, STDP `lambda = 1e-5`**

480 neurons, 80 % of synapses, `corr_F = -0.94`, period CV 0.00034, troughs
deeper than the unreduced network. Reproduced across independent seeds
(-0.943 ± 0.022 and -0.946 ± 0.021) and confirmed at 90 s.

## 1. Synapses: 20 % removable, and that is all

| k (m=0.40) | corr_F | step cost | troughF |
|---|---|---|---|
| 0.250 | -0.946 ± 0.031 | — | 3.69 |
| 0.225 | -0.935 ± 0.031 | +0.011 | 3.41 |
| **0.200** | **-0.946 ± 0.021** | -0.011 | 3.93 |
| 0.175 | -0.909 ± 0.057 | **+0.037** | 4.56 |
| 0.150 | -0.900 ± 0.052 | +0.009 | 5.78 |
| 0.100 | -0.862 ± 0.058 | +0.038 | 5.50 |

k = 0.25, 0.225 and 0.20 are one statistical group. The entire cost is the
single step to 0.175, and the curve flattens again below it. **k = 0.20 is the
floor and there is nothing to gain between 0.20 and 0.25.**

## 2. Neurons: almost nothing is removable

| m (k=0.20) | neurons | corr_F | vs baseline | corr_RG | troughF |
|---|---|---|---|---|---|
| 0.400 | 480 | -0.946 ± 0.021 | +0.007 | -0.945 | 3.93 |
| 0.375 | 450 | -0.936 ± 0.015 | +0.017 | -0.930 | 3.89 |
| 0.350 | 420 | -0.904 ± 0.034 | +0.049 | -0.915 | 4.80 |
| 0.300 | 360 | -0.876 ± 0.053 | +0.077 | -0.874 | 5.73 |

m = 0.375 holds baseline and has the tightest seed SD in the study, but at 90 s
its `corr_RG` reads -0.926 against -0.955 — a gate-independent gap, so the 30
neurons are a real if small loss. Not recommended.

## 3. Why neurons cost more than synapses

The alternation pathway is `rg_e -> in_e -> rg_f`. With `pairwise_bernoulli`
its expected size is `p * N_src * N_tgt` — **linear in k, quadratic in m**.

The alternative explanation, that RG-F simply loses recurrent drive, was tested
and falsified: holding `P_RG_REC * N_RG_F * W` constant by raising
`--w-rg-rec-f` up to 2.5x moved `corr_F` by at most 0.017 in any cell. More
self-excitation raises the DC level, which the gate re-centres away, leaving the
swing unchanged. The limit is structural.

## 4. STDP: no effect on the rhythm, 5x better timing

`k=0.20 / m=0.40`, 90 s, each condition at its own gate.

| | corr_F | corr_RG | period CV | trough_E | trough_F |
|---|---|---|---|---|---|
| off | -0.937 | -0.939 | 0.00173 | 1.81 | 4.29 |
| **on (1e-5)** | -0.936 | -0.938 | **0.00034** | **1.19** | **3.30** |

`corr_F` and `corr_RG` are identical to three decimals — STDP does not change
the rhythm. It cuts period CV **5.1x** and deepens both troughs. Trough depth is
the metric that limits every reduction here, so this is a free gain on the
binding constraint.

At 30 s, STDP appeared to also improve `corr_F` by 0.016. That did not replicate
at 90 s (0.001) and was a transient of unsettled weights. Five of eight 30 s
runs were still drifting; all 90 s runs settled.

## 5. 100 neurons total

8 neurons per population. Both structural levers were swept.

| k | corr_RG (STDP on) | ampF | troughF |
|---|---|---|---|
| 0.25 | -0.240 | 5.22 | 5.48 |
| 0.50 | -0.514 | 5.84 | 4.29 |
| 0.75 | -0.643 | 7.28 | 3.45 |
| **1.00** | **-0.768** | 8.09 | 3.05 |

An earlier attempt at 100 neurons used `k=0.25` and read `corr_RG = -0.29`; that
was a **k artifact, not a size limit**. Raising k to 1.00 lifts alternation by
0.48 — the largest single effect in this study. At small N connectivity buys
more than neurons: 100 neurons at k=1.00 (-0.757) beats 360 neurons at k=0.10
(-0.685).

Inhibition is then exhausted: 3x the weight (`inh_comp` 3.0 → 9.0) buys 0.030 of
`corr_RG` and moves nothing else.

**Ceiling: `corr_RG` ≈ -0.78 at 100 neurons.** k is at 1.00 (all-to-all, no
higher exists) and inhibition is flat.

`corr_F` stalls further back at -0.43, a 0.35 gap that the gate prescription did
not close. RG-F's rate takes only ~10 discrete values at this size (observed:
35.92, 37.21, 38.49, 39.77, 41.06, 42.34, 43.62 Hz — uniformly spaced), so a
sigmoid with `K_F = 30` acts as a comparator rather than interpolating. This is
a readout resolution limit; `grids/n100_kf.tsv` sweeps `K_F` downward to test
whether a shallower gate recovers it.

## 6. Other limits

**Duration.** Frequency drift is under 0.01 % at 60 s. 20 s is enough; 30 s and
90 s agree to within 0.004. Seed variance exceeds any duration effect.

**Speed.** At step period 350 ms the extensor duty cycle rises to 0.67 and its
trough to 3.27 — force cannot decay within the shorter cycle. That is a
muscle-model time constant, not a CPG limit, and no amount of inhibition fixes
it.

## Open

- Every number here is `leg_L`. Inter-leg coordination has never been measured.
- All of this is at 520 ms. The recommended cell has not been re-checked at 780
  or 350 ms.
- `K_F` sweep at 100 neurons is running.
