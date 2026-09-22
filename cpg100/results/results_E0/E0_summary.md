# Summary -- pass B (/Users/jonathan/Documents/CPG-100 neuron/results_E0)

| cell | n | neur | syn | X0_F | corr_F ± SD | Δbase | Δtarget | corr_RG | F−RG | ΔRG ref | corr_F R | troughF | ampF | periodCV | L-R ph | cF t1/t3 | ampF t3/t1 | stable | A1 | A2 | A3 | A4 | A5 | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chk20_ctrl_n480_stdpoff | 8 | 480 | 13100 | 1.256 | -0.872 ± 0.022 | +0.081 |  | -0.634 | -0.238 | +0.304 | -0.846 | 6.70 | 7.73 | 0.00361 | 0.500 | -0.871/-0.881 | 0.98 | 8/8 | -- | -- | ok | -- | ok | fail |
| ctrl_n480_stdpoff | 8 | 480 | 13100 | 1.263 | -0.946 ± 0.026 | +0.007 | -0.006 | -0.949 | +0.004 | -0.011 | -0.892 | 3.65 | 11.79 | 0.00031 | 0.500 | -0.945/-0.947 | 0.99 | 8/8 | ok | ok | ok | ok | ok | PASS |
| ctrl_n480_stdpon | 8 | 480 | 13100 | 1.240 | -0.937 ± 0.023 | +0.016 | +0.003 | -0.946 | +0.009 | +0.000 | -0.881 | 2.88 | 12.73 | 0.00024 | 0.500 | -0.938/-0.953 | 1.08 | 7/8 | ok | ok | ok | ok | ok | PASS |
| base_k100_stdpoff | 8 | 1200 | 118208 | 2.597 | -0.955 ± 0.015 | -0.002 | -0.002 | -0.938 | -0.016 | +0.000 | -0.944 | 0.85 | 15.88 | 0.00004 | 0.500 | -0.953/-0.948 | 1.00 | 6/8 | ok | ok | ok | ok | -- | fail |
| base_k100_stdpon | 8 | 1200 | 118208 | 2.558 | -0.973 ± 0.011 | -0.020 | -0.020 | -0.946 | -0.027 | +0.000 | -0.956 | 0.88 | 15.85 | 0.00001 | 0.500 | -0.972/-0.970 | 1.00 | 6/8 | ok | ok | ok | -- | -- | fail |

A1 |corr_F+0.953|<=0.02   A2 troughF<=4.0   A3 SD(corr_F)<=0.03   A4 corr_RG within 0.02 (|corr_F-corr_RG|)   A5 >= n-1 seeds stable over the run
neur = Izhikevich neurons both legs (built); syn = synapses onto Izhikevich neurons (built);
all metrics leg L except 'corr_F R'; ΔRG ref = corr_RG minus the reference cell's (same STDP);
cF t1/t3 = corr_F in the first/last third of the post-transient window; L-R ph = right-leg extensor onset as a fraction of the left cycle (0.5 = alternation).
