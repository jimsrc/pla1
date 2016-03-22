#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

import numpy as np
import cython_wrapper as cw

a = np.ones(3)

stat = cw.py_run_orbit(a, 0.0, 5.6)
print stat
