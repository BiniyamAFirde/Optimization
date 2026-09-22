#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpg_e2.py -- roadmap node E2: the size ladder 240 -> 180 -> 140 -> 100
======================================================================

Thin layer on cpg_e1.py -> cpg_e0.py -> cpg_small.py v1.1 (nothing copied).
Adds the ladder presets; wiring (--conn-rule indegree, --preserve-input) is set
by the runner. Sizes/leg RGE,RGF,INE,INF,IAINT,MNE,MNF (IaInt built twice):

  n240  24,24,14,19,10,10,9   120/leg   (E1 winner size)
  n180  18,18,11,15, 7, 7,7    90/leg
  n140  14,14, 8,11, 6, 6,5    70/leg
  n100  10,10, 6, 8, 4, 4,4    50/leg   <- the goal: 100 Izhikevich neurons

= the 100-neuron split scaled by N/100 (largest remainder), inhibitory
populations cut last. All: relays 40, k = 0.20, inh-comp 3.75, chunk 100 ms.

With fixed in-degree, K = max(1, round(p * N_src)). At 100 neurons several
pathways hit the floor K = 1 (e.g. RG-E -> InE: p*N = 0.30), i.e. every target
still gets one input where bernoulli would leave ~70 % with none. The run logs
print every projection's in-degree and weight factor -- read them if a rung fails.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cpg_e1  # noqa: E402,F401  (registers n240 in cpg_e0.PRESETS)
import cpg_e0 as W  # noqa: E402

W.E0_VERSION = "cpg_e2 1.0 (cpg_e1 + ladder presets; cpg_small 1.1)"
_COMMON = ["--relay-n", "40", "--k-conn", "0.20", "--inh-comp", "3.75",
           "--inh-comp-f", "1.0", "--chunk-ms", "100"]
for _name, _sizes in (("n180", "18,18,11,15,7,7,7"), ("n140", "14,14,8,11,6,6,5"),
                      ("n100", "10,10,6,8,4,4,4")):
    W.PRESETS[_name] = ["--sizes", _sizes] + _COMMON

if __name__ == "__main__":
    if "--preset" not in sys.argv and "--chunk-ms" not in sys.argv:
        sys.argv[1:1] = ["--chunk-ms", "100"]
    sys.exit(W.main())
