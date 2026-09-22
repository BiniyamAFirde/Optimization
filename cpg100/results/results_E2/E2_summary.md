# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E2)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | A1 | A2 | A3 | A4 | A5 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ctrl_n480_stdpoff | 8 | 480 | 13100 | 1.260 | -0.943 ± 0.023 | +0.010 |  | -0.946 | +0.004 | +0.000 | -0.909 | -0.949 | 1.280 | 3.74 | 11.78 | 0.00034 | 0.500 | -0.946/-0.947 | 0.99 | 7/8 | ok | ok | ok | ok | ok | PASS |
| n240_indeg_pi_stdpoff | 8 | 240 | 5592 | 1.289 | -0.955 ± 0.010 | -0.002 |  | -0.944 | -0.011 | +0.002 | -0.948 | -0.944 | 1.284 | 3.08 | 12.97 | 0.00039 | 0.500 | -0.955/-0.955 | 1.02 | 8/8 | ok | ok | ok | ok | ok | PASS |
| n180_indeg_pi_stdpoff | 8 | 180 | 4072 | 1.263 | -0.900 ± 0.075 | +0.053 |  | -0.931 | +0.031 | +0.015 | -0.936 | -0.932 | 1.295 | 3.74 | 12.26 | 0.00035 | 0.500 | -0.897/-0.903 | 1.03 | 7/8 | -- | ok | -- | -- | ok | fail |
| n180_indeg_stdpoff | 8 | 180 | 4072 | 1.268 | -0.945 ± 0.040 | +0.008 |  | -0.948 | +0.004 | -0.002 | -0.956 | -0.948 | 1.292 | 2.63 | 12.96 | 0.00041 | 0.500 | -0.949/-0.947 | 0.99 | 8/8 | ok | ok | -- | ok | ok | fail |
| n180_indeg_pi_stdpon | 8 | 180 | 4072 | 1.248 | -0.904 ± 0.073 | +0.049 |  | -0.932 | +0.028 | +nan | -0.938 | -0.932 | 1.281 | 3.69 | 12.38 | 0.00030 | 0.500 | -0.896/-0.915 | 1.06 | 8/8 | -- | ok | -- | -- | ok | fail |
| n140_indeg_pi_stdpoff | 8 | 140 | 3108 | 1.258 | -0.902 ± 0.024 | +0.051 |  | -0.912 | +0.009 | +0.035 | -0.911 | -0.909 | 1.266 | 3.14 | 12.61 | 0.00058 | 0.500 | -0.905/-0.902 | 1.05 | 4/8 | -- | ok | ok | ok | ok | fail |
| n140_indeg_stdpoff | 8 | 140 | 3108 | 1.260 | -0.938 ± 0.019 | +0.015 |  | -0.933 | -0.005 | +0.013 | -0.947 | -0.929 | 1.270 | 2.31 | 12.85 | 0.00067 | 0.500 | -0.934/-0.941 | 1.01 | 5/8 | ok | ok | ok | ok | ok | PASS |
| n140_indeg_pi_stdpon | 8 | 140 | 3108 | 1.244 | -0.907 ± 0.023 | +0.046 |  | -0.916 | +0.008 | +nan | -0.915 | -0.909 | 1.254 | 2.99 | 12.92 | 0.00041 | 0.500 | -0.900/-0.918 | 1.03 | 8/8 | -- | ok | ok | ok | ok | fail |
| n100_indeg_pi_stdpoff | 8 | 100 | 2188 | 1.288 | -0.876 ± 0.053 | +0.077 |  | -0.876 | -0.000 | +0.071 | -0.875 | -0.890 | 1.267 | 3.12 | 12.76 | 0.00068 | 0.500 | -0.877/-0.881 | 1.00 | 6/8 | -- | ok | -- | ok | ok | fail |
| n100_indeg_stdpoff | 8 | 100 | 2188 | 1.280 | -0.914 ± 0.046 | +0.039 |  | -0.911 | -0.003 | +0.035 | -0.917 | -0.915 | 1.262 | 2.32 | 12.41 | 0.00101 | 0.500 | -0.909/-0.918 | 1.00 | 7/8 | -- | ok | -- | ok | ok | fail |
| n100_indeg_pi_stdpon | 8 | 100 | 2188 | 1.272 | -0.887 ± 0.040 | +0.066 |  | -0.889 | +0.002 | +nan | -0.880 | -0.891 | 1.259 | 2.85 | 13.03 | 0.00050 | 0.500 | -0.887/-0.900 | 1.02 | 7/8 | -- | ok | -- | ok | ok | fail |

A1 |corr_F+0.953|<=0.02   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_RG within 0.02 (|corr_F-corr_RG|)   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
