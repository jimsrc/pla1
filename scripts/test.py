#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw

AU_in_cm = 1.5e13

ori = np.loadtxt('orientations_isotropic_Nth16_Nph8.in')
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

nB = 10          #--- realizacion de B ---#
npla = 120        # dummy particle id

psim = {
    'rigidity'      : 4.44583E+08,
    'tmax'          : 1e3,
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
    'nHist'         : 100,
    'i'             : npla,
    'j'             : nB,
    'dir_out'       : './',
}


m = cw.mgr()

def set_m(model, nB, pdict, **pother):
    model.set_Bmodel(pdict=pd, nB=nB)
    model.build(**pother)
    #model.runsim(**psim)
    #model.save2file()
    return model

m.set_Bmodel(pdict=pd, nB=nB)
m.build(**pother)

m.runsim(**psim)
m.save2file()

ys = m.ysave

t_sec = m.tsave / m.scl['wc'] # [sec]
pos_x = m.xyz[:,0] * m.scl['rl'] / AU_in_cm

from pylab import plot, show, close, grid, xscale
plot(t_sec, pos_x, '-o')
xscale('log')
grid()
"""
show(); close()
del m
"""
#EOF
