# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E6)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | w drift %/min | A1 | A2 | A3 | A4 | A5 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ctrl_n480_stdpoff | 8 | 480 | 13100 | 1.263 | -0.943 ± 0.024 | +0.010 |  | -0.948 | +0.005 | +0.000 | -0.907 | -0.948 | 1.281 | 3.78 | 11.68 | 0.00031 | 0.500 | -0.942/-0.941 | 1.01 | 8/8 |  | ok | ok | ok | ok | ok | PASS |
| n100_kmin2c_inex15_stdpoff | 8 | 100 | 2353 | 1.189 | -0.965 ± 0.025 | -0.012 |  | -0.948 | -0.017 | -0.000 | -0.952 | -0.946 | 1.184 | 0.99 | 13.86 | 0.00119 | 0.500 | -0.963/-0.965 | 0.99 | 8/8 |  | ok | ok | ok | ok | ok | PASS |
| n100_kmin2c_stdpoff | 8 | 100 | 2353 | 1.247 | -0.949 ± 0.027 | +0.004 |  | -0.941 | -0.008 | +0.007 | -0.929 | -0.937 | 1.238 | 1.32 | 13.37 | 0.00122 | 0.500 | -0.955/-0.945 | 0.98 | 7/8 |  | ok | ok | ok | ok | ok | PASS |
| n100_kmin3c_inex15_stdpoff | 8 | 100 | 2532 | 1.140 | -0.974 ± 0.017 | -0.021 |  | -0.957 | -0.018 | -0.009 | -0.968 | -0.958 | 1.126 | 1.00 | 14.06 | 0.00106 | 0.500 | -0.976/-0.973 | 1.01 | 8/8 |  | -- | ok | ok | ok | ok | fail |
| n100_kmin2c_stdpon | 8 | 100 | 2353 | 1.217 | -0.958 ± 0.020 | +0.015 |  | -0.945 | -0.013 | +nan | -0.938 | -0.941 | 1.212 | 1.17 | 13.59 | 0.00089 | 0.500 | -0.961/-0.960 | 1.01 | 7/8 | 22.73 | ok | ok | ok | ok | ok | PASS |
| n100_kmin3c_inex15_stdpon | 8 | 100 | 2532 | 1.123 | -0.976 ± 0.014 | -0.003 |  | -0.955 | -0.021 | +nan | -0.973 | -0.958 | 1.106 | 0.89 | 14.12 | 0.00081 | 0.500 | -0.973/-0.977 | 1.00 | 8/8 | 22.55 | ok | ok | ok | ok | ok | PASS |
| n100_kmin3c_inex15_stdpon_lam6 | 8 | 100 | 2532 | 1.145 | -0.974 ± 0.015 | -0.001 |  | -0.954 | -0.020 | +nan | -0.969 | -0.958 | 1.124 | 0.91 | 14.05 | 0.00103 | 0.500 | -0.974/-0.972 | 1.00 | 8/8 | 2.81 | ok | ok | ok | ok | ok | PASS |

A1 |corr_F - baseline|<=0.02 (off -0.953, on -0.973)   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_F-corr_RG<=+0.02   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
