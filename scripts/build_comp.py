#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import h5py
from pylab import figure, close
import shared.funcs_post as sfp
from shared.funcs import Bo_parker, Lc_memilia, calc_Rlarmor

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
    for nm in ('xx','yy','zz'):
        mfp = 3.*o['k'+nm]/Lc_adim # [1] normalized units
        lfit[id][nm] = sfp.fit_mfp(tadim, mfp, t_decr=tdecr, **(seeds[nm]))

        ss.update({(id,nm): mfp,})

#--- figs


#EOF
