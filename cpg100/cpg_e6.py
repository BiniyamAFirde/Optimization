#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpg_e6.py -- roadmap node E6: STDP arm + 60 s settle check on the 100-neuron winners
=========================================================================

Thin layer on cpg_e2.py -> cpg_e1.py -> cpg_e0.py -> cpg_small.py v1.2.
No new presets; E3 varies flags the runner passes on top of --preset n100:

  --indeg-min K            minimum inputs per neuron on every Izhikevich-source
                           pathway (cpg_small v1.2). E2 showed all of them sit at
                           K = 1 at 100 neurons, for any k <= 0.6.
  --indeg-min-conserve     split the same total weight over the K inputs
                           (mean drive unchanged, seed variance down)
  --inh-comp G             gain on E -> InE -> F (RG-E->InE and InE->RG-F), recipe 3.75
  --inh-comp-f G           gain on F -> InF -> E (RG-F->InF and InF->RG-E), recipe 1.0
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cpg_e2  # noqa: E402,F401  (registers n240/n180/n140/n100)
import cpg_e0 as W  # noqa: E402

W.E0_VERSION = "cpg_e6 1.0 (cpg_e2 presets; cpg_small 1.2 --indeg-min)"

if __name__ == "__main__":
    if "--preset" not in sys.argv and "--chunk-ms" not in sys.argv:
        sys.argv[1:1] = ["--chunk-ms", "100"]
    sys.exit(W.main())
