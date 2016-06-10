#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import Bparker as bp

x = np.zeros(3, dtype=np.float32)

Rs = np.arange(0.2, 5., 0.05)
x[1] = np.pi/2. # eliptica

Bm = []
for x[0] in Rs:
    B = bp.return_B(x)
    Bm += [ np.square(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]) ]

Bm = np.array(Bm)


