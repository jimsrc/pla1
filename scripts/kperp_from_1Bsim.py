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
dir_src = '../out'
dir_fig = '../figs'
#fname_out = dir_out+'/r.{r:1.2f}_R.{rigidity:1.2e}_eps.{eps_o:1.2e}_NmS.{Nm_slab:04d}_Nm2d.{Nm_2d:04d}_lmin.{lmin:1.1e}.h5'.format(**po)


AUincm = 1.5e13                   # [cm]
po = {}
po['Bo']         = 5e-5                         # [Gauss]
#po['n_modos']    = 128
po['lambda_min'] = 5e-5 #((5e-5)*AUincm)
#--- corregimos input
po['rigidity'] =  1.37352E+08 #4.33306E+07
rl = cw.calc_Rlarmor(po['rigidity'],po['Bo'])   # [cm]
#eps_o = 1.0e-5 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
#po['atol']     = (po['lambda_min']*AUincm)*eps_o/rl
po['rtol']     = 0.0 #1e-6
#fname_inp = './R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.h5'.format(**po)

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

sym = ('o', 's', '^', '*')
#Eps = (3.33e-6, 1e-5, 3.33e-5, 1e-4, 3.33e-4, 1e-3, 3.33e-3, 1e-2,3.33e-2)
eps_o = 1e-5
Ks  = ('kxx', 'kyy', 'kzz')
NMs = [
    # (Nm_slab, Nm_2d, index)
    (256, 256, 0),
    (128, 128, 1),
    (64,  128, 2)
]
lmin = 5e-5
ro   = 0.2

o = {}
for kk in Ks:
    #--- figura
    fig = figure(1, figsize=(6,4))
    ax  = fig.add_subplot(111)
    #for eps_o, i in zip(Eps, range(len(Eps))):
    for Nm_s, Nm_2d, i in NMs:
        po.update({
            'Nm_slab'   : Nm_s,
            'Nm_2d'     : Nm_2d,
            'lmin'      : lmin,
            'eps_o'     : eps_o,
            'rigidity'  : 1.69604E+09,
            'r'         : ro,
        })

        fname_inp = dir_src+'/r.{r:1.2f}_R.{rigidity:1.2e}_eps.{eps_o:1.2e}_NmS.{Nm_slab:04d}_Nm2d.{Nm_2d:04d}_lmin.{lmin:1.1e}.h5'.format(**po)

        o[kk] = ff.get_sqrs(fname_inp)
        tadim = o[kk]['tadim']
        kprof = o[kk][kk]
        #label = '$\epsilon: %2.2e$' % (eps_o)
        label = '$Nm^s, Nm^{2d}: %d, %d$' % (Nm_s, Nm_2d)
        isym =  np.mod(i,len(sym))
        print " i, len(sym), isym: ", i, len(sym), isym;# raw_input()
        msym = sym[isym-1]
        ax.plot(tadim, kprof, '-o', ms=2, lw=0.5, marker=msym, label=label, alpha=0.6, mec='none')

    #ax.set_ylim(1e18, 1e22)#(1e17, 1e21)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('$\Omega t$ [1]')
    ax.set_ylabel('%s [cm2/s]' % kk)
    ax.legend(loc='best', fontsize=6)
    ax.grid()
    #fname_fig = dir_fig+'/{kk}_'.format(kk=kk) +fname_inp.replace(dir_src+'/','').replace('.h5','.png')
    fname_fig = dir_fig+'/{kk}_r.{r:1.2f}_R.{rigidity:1.2e}_eps.{eps_o:1.2e}_lmin.{lmin:1.1e}.png'.format(kk=kk, **po)
    fig.savefig(fname_fig, dpi=300, bbox_inches='tight')
    close(fig)

#EOF
