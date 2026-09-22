#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cpg_e0 as W  # noqa: E402

W.E0_VERSION = "cpg_e1 1.0 (cpg_e0 1.1 + preset n240; cpg_small 1.1)"
W.PRESETS["n240"] = ["--sizes", "24,24,14,19,10,10,9", "--relay-n", "40",
                     "--k-conn", "0.20", "--inh-comp", "3.75", "--inh-comp-f", "1.0",
                     "--chunk-ms", "100"]

if __name__ == "__main__":
    if "--preset" not in sys.argv and "--chunk-ms" not in sys.argv:
        sys.argv[1:1] = ["--chunk-ms", "100"]
    sys.exit(W.main())
