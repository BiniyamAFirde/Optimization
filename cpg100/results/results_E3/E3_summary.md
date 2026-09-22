# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E3)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | A1 | A2 | A3 | A4 | A5 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ctrl_n480_stdpoff | 8 | 480 | 13100 | 1.260 | -0.943 ± 0.023 | +0.010 |  | -0.946 | +0.004 | +0.000 | -0.909 | -0.949 | 1.280 | 3.74 | 11.78 | 0.00034 | 0.500 | -0.946/-0.947 | 0.99 | 7/8 | ok | ok | ok | ok | ok | PASS |
| n140_indeg_stdpoff | 8 | 140 | 3108 | 1.260 | -0.938 ± 0.019 | +0.015 |  | -0.933 | -0.005 | +0.013 | -0.947 | -0.929 | 1.270 | 2.31 | 12.85 | 0.00067 | 0.500 | -0.934/-0.941 | 1.01 | 5/8 | ok | ok | ok | ok | ok | PASS |
| n100_indeg_stdpoff | 8 | 100 | 2188 | 1.280 | -0.914 ± 0.046 | +0.039 |  | -0.911 | -0.003 | +0.035 | -0.917 | -0.915 | 1.262 | 2.32 | 12.41 | 0.00101 | 0.500 | -0.909/-0.918 | 1.00 | 7/8 | -- | ok | -- | ok | ok | fail |
| n100_kmin2c_stdpoff | 8 | 100 | 2353 | 1.238 | -0.952 ± 0.028 | +0.001 |  | -0.941 | -0.010 | +0.005 | -0.930 | -0.939 | 1.234 | 1.64 | 13.14 | 0.00116 | 0.500 | -0.952/-0.952 | 1.01 | 7/8 | ok | ok | ok | ok | ok | PASS |
| n100_kmin3_stdpoff | 8 | 100 | 2532 | 1.038 | -0.964 ± 0.019 | -0.011 |  | -0.936 | -0.028 | +0.010 | -0.964 | -0.940 | 1.041 | 0.62 | 14.61 | 0.00123 | 0.500 | -0.962/-0.973 | 0.99 | 7/8 | ok | ok | ok | -- | ok | fail |
| n100_kmin3c_inex05_stdpoff | 8 | 100 | 2532 | 1.296 | -0.902 ± 0.063 | +0.051 |  | -0.901 | -0.001 | +0.045 | -0.881 | -0.907 | 1.257 | 2.51 | 11.69 | 0.00110 | 0.500 | -0.904/-0.909 | 1.03 | 6/8 | -- | ok | -- | ok | ok | fail |
| n100_kmin3c_inex15_stdpoff | 8 | 100 | 2532 | 1.149 | -0.970 ± 0.020 | -0.017 |  | -0.956 | -0.014 | -0.009 | -0.975 | -0.958 | 1.136 | 0.94 | 13.98 | 0.00101 | 0.500 | -0.973/-0.972 | 1.01 | 7/8 | ok | ok | ok | ok | ok | PASS |
| n100_kmin3c_infx2_stdpoff | 8 | 100 | 2532 | 1.228 | -0.950 ± 0.040 | +0.003 |  | -0.947 | -0.003 | -0.001 | -0.930 | -0.952 | 1.181 | 1.27 | 13.50 | 0.00105 | 0.500 | -0.947/-0.952 | 1.01 | 8/8 | ok | ok | -- | ok | ok | fail |
| n100_kmin3c_stdpoff | 8 | 100 | 2532 | 1.218 | -0.952 ± 0.036 | +0.001 |  | -0.949 | -0.003 | -0.003 | -0.934 | -0.949 | 1.184 | 1.45 | 13.35 | 0.00107 | 0.500 | -0.955/-0.948 | 1.00 | 6/8 | ok | ok | -- | ok | ok | fail |
| n100_kmin3c_stdpon | 8 | 100 | 2532 | 1.211 | -0.957 ± 0.031 | -0.004 |  | -0.949 | -0.008 | +nan | -0.938 | -0.952 | 1.175 | 1.20 | 13.57 | 0.00085 | 0.500 | -0.955/-0.963 | 1.03 | 8/8 | ok | ok | -- | ok | ok | fail |

A1 |corr_F+0.953|<=0.02   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_RG within 0.02 (|corr_F-corr_RG|)   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
