# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E8)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | w drift %/min | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n100_kmin2c_g350_stdpoff | 8 | 100 | 2353 | 1.246 | -0.892 ± 0.023 | +0.026 |  | -0.909 | +0.016 | -0.013 | -0.864 | -0.895 | 1.246 | 2.55 | 11.20 | 0.00185 | 0.500 | -0.892/-0.898 | 0.99 | 7/8 |  | -- | ok | ok | ok | ok | -- | ok | ok | fail |
| n100_kmin2c_g520_stdpoff | 8 | 100 | 2353 | 1.247 | -0.949 ± 0.027 | +0.006 |  | -0.941 | -0.008 | -0.046 | -0.929 | -0.937 | 1.238 | 1.32 | 13.37 | 0.00122 | 0.500 | -0.955/-0.945 | 0.98 | 7/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_kmin2c_g780_stdpoff | 8 | 100 | 2353 | 1.240 | -0.963 ± 0.017 | +0.031 |  | -0.830 | -0.133 | +0.065 | -0.943 | -0.822 | 1.229 | 1.56 | 13.28 | 0.00157 | 0.500 | -0.962/-0.964 | 1.01 | 8/8 |  | -- | ok | ok | ok | ok | -- | ok | -- | fail |
| n100_kmin3c_inex15_lam6_g350_stdpon | 8 | 100 | 2532 | 1.164 | -0.905 ± 0.038 | +0.014 |  | -0.918 | +0.013 | +nan | -0.911 | -0.922 | 1.147 | 2.01 | 11.48 | 0.00181 | 0.500 | -0.905/-0.901 | 1.00 | 7/8 | 2.72 | ok | ok | -- | ok | ok | ok | ok | ok | fail |
| n100_kmin3c_inex15_lam6_g520_stdpon | 8 | 100 | 2532 | 1.145 | -0.974 ± 0.015 | -0.020 |  | -0.954 | -0.020 | +nan | -0.969 | -0.958 | 1.124 | 0.91 | 14.05 | 0.00103 | 0.500 | -0.974/-0.972 | 1.00 | 8/8 | 2.81 | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| n100_kmin3c_inex15_lam6_g780_stdpon | 8 | 100 | 2532 | 1.120 | -0.978 ± 0.017 | +0.016 |  | -0.873 | -0.105 | +nan | -0.974 | -0.885 | 1.106 | 0.96 | 14.38 | 0.00177 | 0.500 | -0.978/-0.981 | 1.00 | 8/8 | 2.89 | ok | ok | ok | ok | ok | ok | ok | -- | fail |
| base_k100_g350_stdpoff | 8 | 1200 | 118208 | 2.895 | -0.918 ± 0.014 | +0.000 |  | -0.895 | -0.023 | +0.000 | -0.917 | -0.894 | 2.952 | 2.52 | 13.90 | 0.00202 | 0.500 | -0.913/-0.929 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g520_stdpoff | 8 | 1200 | 118208 | 2.597 | -0.955 ± 0.013 | +0.000 |  | -0.938 | -0.017 | -0.043 | -0.942 | -0.936 | 2.601 | 0.83 | 15.89 | 0.00000 | 0.500 | -0.951/-0.962 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |
| base_k100_g780_stdpoff | 8 | 1200 | 118208 | 2.260 | -0.994 ± 0.003 | +0.000 |  | -0.930 | -0.064 | -0.035 | -0.993 | -0.923 | 2.257 | 0.25 | 17.06 | 0.00194 | 0.500 | -0.992/-0.995 | 1.00 | 8/8 |  | ok | ok | ok | ok | ok | ok | ok | ok | PASS |

A1 |corr_F - baseline|<=0.02 (off -0.953, on -0.973)   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_F-corr_RG<=+0.02   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
