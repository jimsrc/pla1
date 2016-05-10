#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from pylab import (
    plot, scatter, show, close, grid, figure, xticks, 
    xscale, yscale, xlabel, ylabel, bar, find, hist, 
    contourf, savefig
)
import numpy as np
from os.path import isfile, isdir
import src.cython_wrapper as cw
import shared.funcs as ff

"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01
1e10    1.0898E+10   4.853544E-02   9.963142E-01
"""
AUincm = 1.5e13                   # [cm]
po = {}
po['Bo']         = 5e-5                         # [Gauss]
po['n_modos']    = 128
po['lambda_min'] = 5e-5 #((5e-5)*AUincm)
#--- corregimos input
po['rigidity'] =  1.37352E+08 #4.33306E+07
rl = cw.calc_Rlarmor(po['rigidity'],po['Bo'])   # [cm]
eps_o = 1.0e-5 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
po['atol']     = (po['lambda_min']*AUincm)*eps_o/rl
po['rtol']     = 0.0 #1e-6
fname_inp = './R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.h5'.format(**po)

#print fname_inp
"""
o = ff.get_sqrs(fname_inp)
nplas = o['nplas']
ntime = o['ntime']
kxx = o['kxx']
kyy = o['kyy']
kzz = o['kzz']
tadim = o['tadim']
tdim = o['tdim']
wc = o['wc']
rl = o['rl']
"""


Eps = (3.33e-6, 1.0e-5, 3.33e-5, 1.0e-4)
#Eps = (1.0e-5, 3.33e-5, 1.0e-4)
Ks  = ('kxx', 'kyy', 'kzz')
o = {}
for kk in Ks:
    #--- figura
    fig = figure(1, figsize=(6,4))
    ax  = fig.add_subplot(111)
    for eps_o in Eps:
        po['atol']     = (po['lambda_min']*AUincm)*eps_o/rl
        fname_inp = './R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.h5'.format(**po)
        o[kk] = ff.get_sqrs(fname_inp)
        tadim = o[kk]['tadim']
        kprof = o[kk][kk]
        label = '$\epsilon: %2.2e, atol:%2.1e$' % (eps_o, po['atol'])
        ax.plot(tadim, kprof, '-o', ms=3, label=label)

    ax.set_ylim(1e17, 1e21)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('%s [cm2/s]' % kk)
    ax.legend(loc='best')
    ax.grid()
    fname_fig = './R.{rigidity:1.2e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.png'.format(**po)
    fig.savefig(fname_fig, dpi=135)
    close(fig)

#EOF
