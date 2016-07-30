#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
#import shared.funcs as ff
from pylab import pause
#--- parameters
from params import (
    nB, pd, psim, pother, mu, ph, AU_in_cm
)
#from mpi4py import MPI
from h5py import File as h5
from os.path import isfile, isdir
import os, sys
from glob import glob
#from Bparker.Bparker import return_B as Bparker

"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01  (*)
1e10    1.0898E+10   4.853544E-02   9.963142E-01
--- Bparker model
(Ek=1e9eV)
r[AU]    B[nT]       Rl[AU]         Lc[AU]      Rl/Lc   Rl/(5e-5AU)
0.2      91.0542857  4.147812E-04   0.0044549   0.09    8.30
0.3      41.4298571  9.116035E-04   0.0053034   0.17    18.23
0.4      24.041      1.570966E-03   0.0060017   0.26    31.42
0.5      15.972      2.364613E-03   0.0066061   0.36    47.29
0.7      8.897       4.244982E-03   0.0076345   0.56    84.90
0.8      7.14642857  5.284822E-03   0.0080857   0.65    105.70
1.0      5.0         7.553521E-03   0.0089      0.85    151.07
2.0      1.99653571  1.891657E-02   0.0119904   1.58    378.33
"""
#--- set B-turbulence model
Lc_slab = 0.01 # ... (maria emilia) # [AU]
Rl = 1. # [AU] #cw.calc_Rlarmor(psim['rigidity'],pd['Bo']) #[cm]
pd.update({
    'Nm_slab'       : 64,
    'Nm_2d'         : 64,
    'lmin_s'        : 5e-5/Rl, #[lmin_s/Lc_slab] #(5e-5)*AU_in_cm,
    'lmax_s'        : 1.0/Rl,  #[lmax_s/Lc_slab] #(1.0) *AU_in_cm,
    'lmin_2d'       : 5e-5/Rl, #[lmin_2d/Lc_slab] #(5e-5)*AU_in_cm,
    'lmax_2d'       : 1.0/Rl,  #[lmax_2d/Lc_slab] #(1.0) *AU_in_cm,
    'Lc_slab'       : Lc_slab/Rl,  # in units of Larmor-radii
    'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
    'sigma_Bo_ratio': 0.3, # [1] fluctuation energy
    #'Bo'            : Bo, # [Gauss]
    'ratio_slab'    : 0.2, # [1] (energy_slab)/(energy_total)
})

#--- corregimos input
#Rl_o_LcSlab = ... #(ratio of Larmor radii and Lc_slab)
#eps_o = 3.33e-5 #3.33e-6 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
#lmin = np.min([pd['lmin_s'], pd['lmin_2d']]) # [cm] smallest turb scale

#--- call simulator
m = cw.mgr()
m.set_Bmodel(pdict=pd, nB=nB)
print m.Bxyz([0.,0.,0.])


#EOF
