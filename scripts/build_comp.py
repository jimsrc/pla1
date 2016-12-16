#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import h5py, argparse
from pylab import figure, close
import shared.funcs_post as sfp
from shared.funcs import Bo_parker, Lc_memilia, calc_Rlarmor

# args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-fig', '--figname', 
type=str,
help='figname.',
)
pa = parser.parse_args()

ro = 1 # [AU]
#--- seeds for fit
seeds = {
'xx': {'seed_b': 0.1, 'seed_m': 190.},
'yy': {'seed_b': 0.1, 'seed_m': 190.},
'zz': {'seed_b': 30., 'seed_m': -1.6e4},
}
tdecr = 5e2
#fname_inp = '../out/auger_001.h5'
fs = {
1: ['../out/auger_001.h5', 1.37352E+08],
2: ['../out/auger_002.h5', 1.69604E+09],
}
lfit = {}
ss   = {}

for id in fs.keys():
    fname_inp = fs[id][0]
    o = sfp.get_sqrs(fname_inp)
    f = h5py.File(fname_inp, 'r')
    rig = fs[id][1] #1.37352E+08 # [V] Ek=1e7eV
    Rl, v, omg = calc_Rlarmor(
            rigidity=rig, 
            Bo=Bo_parker(r=ro),
            full=True
            ) # [cm], [cm/s], [s^-1]
    Lc = Lc_memilia(r=ro)

    tadim = o['tadim'] # [1]
    lfit[id] = {}
    Lc_adim = f['psim/Lc_slab'].value # [1]
    ss.update({(id,'t'): tadim/omg}) # [seg]
    for nm in ('xx','yy','zz'):
        mfp = 3.*o['k'+nm]/Lc_adim # [1] normalized units
        lfit[id][nm] = sfp.fit_mfp(tadim, mfp, t_decr=tdecr, **(seeds[nm]))
        ss.update({(id,nm): mfp,})

#--- figs
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

#id=1
for id in (1,2):
    lperp = .5*(ss[(id,'xx')] + ss[(id,'yy')]) # [1]
    time  = ss[(id,'t')]  # [sec]
    lperp_fit = 0.5*(lfit[id]['xx'] + lfit[id]['yy'])
    ax.plot(time/3600., lperp/lperp_fit, '-o', ms=3, label='$R=%2.2e GV$'%(fs[id][1]/1e9))

ax.set_xlabel('time [hr]')
ax.set_ylabel('$\lambda_\perp / \lambda_{\perp O}$ [1]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.legend(loc='best')

fig.savefig(pa.figname, dpi=135, bbox_inches='tight')
close(fig)

#EOF
