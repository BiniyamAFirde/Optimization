# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E1)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | corr_RG R | X0_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | A1 | A2 | A3 | A4 | A5 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ctrl_n480_stdpoff | 8 | 480 | 13100 | 1.260 | -0.943 ± 0.023 | +0.010 |  | -0.946 | +0.004 | +0.000 | -0.909 | -0.949 | 1.280 | 3.74 | 11.78 | 0.00034 | 0.500 | -0.946/-0.947 | 0.99 | 7/8 | ok | ok | ok | ok | ok | PASS |
| n240_bern_pi_stdpoff | 8 | 240 | 5509 | 1.330 | -0.839 ± 0.113 | +0.114 |  | -0.898 | +0.059 | +0.048 | -0.868 | -0.903 | 1.287 | 3.91 | 9.60 | 0.00043 | 0.500 | -0.828/-0.855 | 1.03 | 8/8 | -- | ok | -- | -- | ok | fail |
| n240_bern_stdpoff | 8 | 240 | 5509 | 1.337 | -0.835 ± 0.109 | +0.118 |  | -0.890 | +0.056 | +0.056 | -0.869 | -0.894 | 1.301 | 4.28 | 9.01 | 0.00048 | 0.500 | -0.821/-0.840 | 0.96 | 5/8 | -- | -- | -- | -- | ok | fail |
| n240_indeg_pi_stdpoff | 8 | 240 | 5592 | 1.289 | -0.955 ± 0.010 | -0.002 |  | -0.944 | -0.011 | +0.002 | -0.948 | -0.944 | 1.284 | 3.08 | 12.97 | 0.00039 | 0.500 | -0.955/-0.955 | 1.02 | 8/8 | ok | ok | ok | ok | ok | PASS |
| n240_indeg_stdpoff | 8 | 240 | 5592 | 1.289 | -0.972 ± 0.006 | -0.019 |  | -0.955 | -0.018 | -0.008 | -0.959 | -0.953 | 1.289 | 2.32 | 13.64 | 0.00040 | 0.500 | -0.978/-0.969 | 1.00 | 7/8 | ok | ok | ok | ok | ok | PASS |
| n240_indeg_pi_stdpon | 8 | 240 | 5592 | 1.273 | -0.952 ± 0.011 | +0.001 |  | -0.942 | -0.011 | +nan | -0.938 | -0.942 | 1.267 | 2.88 | 13.30 | 0.00027 | 0.500 | -0.948/-0.962 | 1.04 | 8/8 | ok | ok | ok | ok | ok | PASS |

A1 |corr_F+0.953|<=0.02   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_RG within 0.02 (|corr_F-corr_RG|)   A5 cell mean: dcorr_F(t3-t1)<=+0.02, ampF t3/t1>=0.90, cycles ok in >=n-1 seeds ('stable' = seeds passing all three alone)
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
