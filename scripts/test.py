#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
import shared.funcs as ff

#-----------------------------------------
AU_in_cm = 1.5e13

#ori = np.loadtxt('orientations_isotropic_Nth16_Nph8.in')
ori = np.loadtxt('orientations_122.in')
mu  = np.cos(ori[:,0])
ph  = ori[:,1]

#--- B-field parameters only
pd = {
    # parametros fisicos
    'n_modos'       : 128,
    'lambda_min'    : ((5e-5)*AU_in_cm),
    'lambda_max'    : 1.0*AU_in_cm,
    'Lc_slab'       : 0.01*AU_in_cm,
    'Lc_2d'         : 0.01*AU_in_cm,
    'sigma_Bo_ratio': 1.0,
    'percent_slab'  : 0.2,
    'percent_2d'    : 0.8,
    'Bo'            : 5e-5,   # [Gauss]
    # semillas
    'sem_slab0'     : 17,
    'sem_slab1'     : 101,
    'sem_slab2'     : 33,
    'sem_two0'      : 14,
    'sem_two1'      : 79,
}

nB = 0          #--- realizacion de B ---#
npla = 0        # dummy particle id

"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01
1e10    1.0898E+10   4.853544E-02   9.963142E-01
"""
psim = {
    'rigidity'      : 4.33306E+07,
    'tmax'          : 1e4, #4e4,
    'FracGyroperiod': 5e-2,
    'hmin'          : 0.0,
    'mu'            : mu[npla],
    'ph'            : ph[npla],
    'rtol'          : 1e-5,
    'atol'          : 1e-5,
}

pother = {
    'str_timescale' : 'mixed',
    'nsave'         : 20,
    'tmaxHistTau'   : 150,
    'nHist'         : 150,
    'i'             : npla,
    'j'             : nB,
    'dir_out'       : './',
}


m = cw.mgr()

m.set_Bmodel(pdict=pd, nB=nB)
m.build(**pother)

m.SetSim(**psim)
m.RunSim()

fname_out = 'test.h5'
ff.save_to_h5(m, fname_out)

h_step = m.HistStep
h_seq  = m.HistSeq
#minst = m.MinStep

#m.runsim(**psim)
#m.set_sim(**psim)
#print " ---->>> ok.... "
#m.runsim_partial(4e2)
"""
#m.set_all(pdict=pd, nB=nB, pother=pother)

m.runsim(**psim)
#m.save2file()

ys = m.ysave

t_sec = m.tsave / m.scl['wc'] # [sec]
pos_x = m.xyz[:,0] * m.scl['rl'] / AU_in_cm

scl   = m.scl
wc    = scl['wc']
err   = m.err
tadim = t_sec*wc

import numpy as np
from h5py import File as h5
from pylab import plot, show, close, grid, xscale, yscale, xlabel, ylabel

fname_out = '_prueba.h5'
f = h5(fname_out, 'w')
f['t_sec'] = t_sec
f['pos_x'] = pos_x
f['ysave'] = ys
f['wc'] = wc
f['err'] = err
f['tadim'] = tadim
f.close()
"""

"""
# abremos output del repo de massive-orbits:
from h5py import File as h5
f2 = h5('/home/jim/simulacion/pla_stochastic/composite_turbulence/massive_orbits/output/Ek.1.0e+06eV_rtol.1e-06/out.h5')
x_f2 = f2['B%02d/pla%03d/xyz'%(nB,npla)][0,:]

# ploteamos para comparar
plot(tadim, x_f2,  '-o')
plot(tadim, pos_x, '-o')
grid()
xscale('log')

show(); close()
"""

"""
from pylab import plot, show, close, grid, xscale
plot(t_sec, pos_x, '-o')
xscale('log')
grid()

show(); close()
del m
"""
#EOF
