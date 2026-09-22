# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E8g)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | w drift %/min | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n100_final10_g350_stdpoff | 8 | 100 | 2710 | 1.218 | -0.916 ± 0.025 | +0.002 |  | -0.914 | -0.003 | -0.019 | -0.903 | -0.913 | 1.241 | 2.38 | 11.33 | 0.00172 | 0.499 | -0.917/-0.913 | 1.01 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final10_g520_stdpoff | 8 | 100 | 2710 | 1.210 | -0.979 ± 0.008 | -0.024 |  | -0.953 | -0.026 | -0.058 | -0.960 | -0.950 | 1.230 | 1.12 | 13.69 | 0.00105 | 0.500 | -0.980/-0.982 | 1.00 | 8/8 |  | -- | ok | ok | ok | ok | ok | ok | ok | fail |
| n100_final10_g780_stdpoff | 8 | 100 | 2710 | 1.881 | -0.991 ± 0.003 | +0.002 |  | -0.900 | -0.091 | -0.005 | -0.990 | -0.900 | 1.910 | 0.38 | 15.79 | 0.00151 | 0.500 | -0.992/-0.991 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | -- | fail |
| n100_final10_g350_stdpon | 8 | 100 | 2710 | 1.216 | -0.921 ± 0.024 | -0.002 |  | -0.916 | -0.005 | +nan | -0.906 | -0.912 | 1.241 | 2.39 | 11.27 | 0.00164 | 0.499 | -0.923/-0.920 | 1.00 | 5/8 | 2.72 | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final10_g520_stdpon | 8 | 100 | 2710 | 1.210 | -0.980 ± 0.006 | -0.025 |  | -0.953 | -0.027 | +nan | -0.966 | -0.951 | 1.233 | 1.07 | 13.70 | 0.00099 | 0.500 | -0.981/-0.980 | 1.01 | 8/8 | 2.80 | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_final10_g780_stdpon | 8 | 100 | 2710 | 1.868 | -0.990 ± 0.002 | +0.003 |  | -0.897 | -0.094 | +nan | -0.989 | -0.899 | 1.906 | 0.39 | 15.77 | 0.00156 | 0.500 | -0.989/-0.991 | 1.00 | 8/8 | 2.86 | ok | ok | ok | ok | ok | ok | ok | -- | fail |
| base_k100_g350_stdpoff | 8 | 1200 | 118208 | 2.895 | -0.918 ± 0.014 | +0.000 |  | -0.895 | -0.023 | +0.000 | -0.917 | -0.894 | 2.952 | 2.52 | 13.90 | 0.00202 | 0.500 | -0.913/-0.929 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g520_stdpoff | 8 | 1200 | 118208 | 2.597 | -0.955 ± 0.013 | +0.000 |  | -0.938 | -0.017 | -0.043 | -0.942 | -0.936 | 2.601 | 0.83 | 15.89 | 0.00000 | 0.500 | -0.951/-0.962 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g780_stdpoff | 8 | 1200 | 118208 | 2.260 | -0.994 ± 0.003 | +0.000 |  | -0.930 | -0.064 | -0.035 | -0.993 | -0.923 | 2.257 | 0.25 | 17.06 | 0.00194 | 0.500 | -0.992/-0.995 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |

A1 |corr_F - baseline|<=0.02 (off -0.953, on -0.973)   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_F-corr_RG<=+0.02   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
