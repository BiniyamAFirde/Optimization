# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E8f)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | w drift %/min | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n100_final_g350_stdpoff | 8 | 100 | 2710 | 1.146 | -0.932 ± 0.015 | -0.013 |  | -0.921 | -0.011 | -0.026 | -0.925 | -0.922 | 1.180 | 2.09 | 11.46 | 0.00165 | 0.499 | -0.928/-0.931 | 0.99 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final_g520_stdpoff | 8 | 100 | 2710 | 1.132 | -0.987 ± 0.007 | -0.033 |  | -0.955 | -0.032 | -0.060 | -0.981 | -0.957 | 1.155 | 0.87 | 13.99 | 0.00096 | 0.500 | -0.986/-0.987 | 1.00 | 8/8 |  | -- | ok | ok | ok | ok | -- | ok | ok | fail |
| n100_final_g780_stdpoff | 8 | 100 | 2710 | 1.713 | -0.993 ± 0.002 | +0.001 |  | -0.919 | -0.073 | -0.024 | -0.992 | -0.924 | 1.752 | 0.27 | 15.93 | 0.00156 | 0.500 | -0.993/-0.993 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final_g350_stdpon | 8 | 100 | 2710 | 1.139 | -0.927 ± 0.019 | -0.009 |  | -0.921 | -0.007 | +nan | -0.922 | -0.921 | 1.182 | 2.10 | 11.52 | 0.00161 | 0.499 | -0.925/-0.932 | 1.00 | 7/8 | 2.72 | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final_g520_stdpon | 8 | 100 | 2710 | 1.129 | -0.987 ± 0.005 | -0.032 |  | -0.958 | -0.029 | +nan | -0.982 | -0.957 | 1.159 | 0.87 | 14.00 | 0.00094 | 0.500 | -0.988/-0.987 | 1.00 | 8/8 | 2.81 | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final_g780_stdpon | 8 | 100 | 2710 | 1.693 | -0.992 ± 0.003 | +0.002 |  | -0.921 | -0.071 | +nan | -0.993 | -0.923 | 1.747 | 0.28 | 15.92 | 0.00155 | 0.500 | -0.991/-0.992 | 1.00 | 8/8 | 2.86 | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g350_stdpoff | 8 | 1200 | 118208 | 2.895 | -0.918 ± 0.014 | +0.000 |  | -0.895 | -0.023 | +0.000 | -0.917 | -0.894 | 2.952 | 2.52 | 13.90 | 0.00202 | 0.500 | -0.913/-0.929 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g520_stdpoff | 8 | 1200 | 118208 | 2.597 | -0.955 ± 0.013 | +0.000 |  | -0.938 | -0.017 | -0.043 | -0.942 | -0.936 | 2.601 | 0.83 | 15.89 | 0.00000 | 0.500 | -0.951/-0.962 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g780_stdpoff | 8 | 1200 | 118208 | 2.260 | -0.994 ± 0.003 | +0.000 |  | -0.930 | -0.064 | -0.035 | -0.993 | -0.923 | 2.257 | 0.25 | 17.06 | 0.00194 | 0.500 | -0.992/-0.995 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |

A1 |corr_F - baseline|<=0.02 (off -0.953, on -0.973)   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_F-corr_RG<=+0.02   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
