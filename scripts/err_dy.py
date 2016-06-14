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
from h5py import File as h5

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

AUincm = 1.5e13                   # [cm]
po = {}
#po['Lc_slab']    = 0.01*AUincm
po['Bo']         = 5e-5                         # [Gauss]
po['n_modos']    = 128
po['lambda_min'] = 5e-5 #((5e-5)*AUincm)
#--- corregimos input
po['rigidity'] =  1.37352E+08 #4.33306E+07
rl = cw.calc_Rlarmor(po['rigidity'],po['Bo'])   # [cm]
#eps_o = 1.0e-5 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)

sym = ('o', 's', '^', '*')
#Eps = (3.33e-6, 1e-5, 3.33e-5, 1e-4, 3.33e-4, 1e-3, 3.33e-3, 1e-2,3.33e-2)
eps_o = 3.33e-4
NMs = [
    # (Nm_slab, Nm_2d, index)
    (256, 256, 0),
    (128, 128, 1),
    (64,  128, 2)
]
lmin = 5e-5
ro   = 0.2

#neps = len(Eps)
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)
ax2 = ax.twiny()
#for eps_o, ie in zip(Eps, range(neps)):
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
    #fname_inp = './R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.h5'.format(**po)
    #for eps_o, i in zip(Eps, range(len(Eps))):
    f = h5(fname_inp, 'r')

    wc = f['pla000/scl_wc'].value             # [s-1]
    vp = f['pla000/scl_vel'].value            # [cm/s]
    rl = f['pla000/scl_rl'].value             # [cm]
    PNAMES = f.keys()
    Np     = len(PNAMES)

    hmg = ff.Hmgr(f, nbin=1000)
    hmg.get_hstep_extremes()

    for pnm, ip in zip(PNAMES, range(Np)):
        if not pnm.startswith('pla'):
            continue

        h  = f[pnm+'/HistStep/HStep'].value
        hx = f[pnm+'/HistStep/bins_StepPart'].value
        hmg.pile_to_hist(hx, h)

    isym = np.mod(i,len(sym))
    msym = sym[isym-1]
    opt = {'ms': 3, 'mec':'none', 'marker': msym,'ls':''}
    #label = '$\epsilon: %1.1e$' % eps_o #+\
    label = '$Nm^s, Nm^{2d}: %d, %d$' % (Nm_s, Nm_2d)
    dRbin = (hmg.hbin/wc)*vp/(po['lambda_min']*AUincm)
    ax2.plot(dRbin, hmg.h, label=label, **opt)
    ax2.set_xlim(dRbin[0], dRbin[-1])
    ax2.set_xlabel('$\Delta r/\lambda_{min}$')
    ax2.set_xscale('log')
    ax.plot(hmg.hbin, hmg.h, label=label, **opt)
    ax.grid()
    ax.legend(loc='best', fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$\Omega dt$')
    ax.set_ylabel('#')
    ax.set_xlim(hmg.hbin[0], hmg.hbin[-1])
    #ax.set_ylim(1.0, 1e6)
    #print " ---> generating: " + fname_fig
    f.close()

#fname_fig = './R.{rigidity:1.2e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.png'.format(**po)
fname_fig = dir_fig+'/errdy_r.{r:1.2f}_R.{rigidity:1.2e}_eps.{eps_o:1.2e}_lmin.{lmin:1.1e}.png'.format(**po)
fig.savefig(fname_fig, dpi=200, bbox_inches='tight')
close(fig)

#EOF
