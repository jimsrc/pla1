#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
import funcs as ff
import os
from os.path import isfile, isdir
#--- parameters
from params import (
    nB, pd, psim, pother, mu, ph, AU_in_cm
)

"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01
1e10    1.0898E+10   4.853544E-02   9.963142E-01
"""

nB = 0          #--- realizacion de B ---#
npla = 0        # dummy particle id
#--- set B-turbulence model
pd['n_modos']    = 128
pd['lambda_min'] = ((5e-5)*AU_in_cm)
#--- corregimos input
psim['rigidity'] = 1.69604E+09
rl = cw.calc_Rlarmor(psim['rigidity'],pd['Bo']) #[cm]
eps_o = 4e-5
psim['atol']     = pd['lambda_min']*eps_o/rl
psim['rtol']     = 0.0 #1e-6
psim['tmax']     = 4e4 #0.3e4 #4e4
print " ----> simulation parameters:\n", psim

m = cw.mgr()

m.set_Bmodel(pdict=pd, nB=nB)
m.build(**pother)

m.SetSim(**psim)
m.RunSim()

fname_out = 'test.h5'
ff.save_to_h5(m, fname_out)

h_step = m.HistStep
h_seq  = m.HistSeq

#--- output
po = {}
po.update(psim)
po.update(pd)
po['lambda_min'] /= AU_in_cm
dir_out = './err_out'
fname_out = dir_out+'/R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.h5'.format(**po)

#--- save to file
f = ff.save_to_h5(m, fname_out, file=None, close=False)
f['HistStep/htot_bin'] = h_step[:,0]
f['HistStep/htot_cts'] = h_step[:,1]
f['HistStep/hrel_bin'] = h_step[:,2]
f['HistStep/hrel_cts'] = h_step[:,3]
f['HistSeq/seq_bin']   = h_seq[:,0]
f['HistSeq/seq_cts']   = h_seq[:,1]
f.close()

#m.runsim(**psim)
#m.set_sim(**psim)
#print " ---->>> ok.... "
#m.runsim_partial(4e2)
#EOF
