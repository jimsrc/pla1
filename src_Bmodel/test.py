#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
import Bmodel
#from Bparker.Bparker import return_B as Bparker_vector
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia
#--- parameters
from params import (
    nB, pd, psim, pother, mu, ph, AU_in_cm
)


AUincm = 1.5e13     # [cm]
xyz = np.zeros(3, dtype=np.float32)

m  = Bmodel.Bmodel() # model slab/2D
ro = 1.0                   # [AU]
Lc_slab = Lc_memilia(r=ro) # [AU]
Bo      = Bo_parker(r=ro)  # [G]
Rl = calc_Rlarmor(
    rigidity=1e9,       # [V] rigidity
    Bo=Bo,              # [G] Bo-field
    )/AUincm            # [AU]
#--- set B-turbulence model
pd.update({
'Nm_slab'       : 128,
'Nm_2d'         : 128,
'lmin_s'        : 5e-5/Rl, #[lmin_s/Rl] 
'lmax_s'        : 1.0/Rl,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-5/Rl, #[lmin_2d/Rl] 
'lmax_2d'       : 1.0/Rl,  #[lmax_2d/Rl] 
'Lc_slab'       : Lc_slab/Rl,  # in units of Larmor-radii
'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
'sigma_Bo_ratio': 0.3, # [1] fluctuation energy
'ratio_slab'    : 0.2, # [1] (energy_slab)/(energy_total)
})
m.set_Bmodel(pdict=pd, nB=nB)
print m.Bxyz(xyz)


"""
Rs = np.arange(0.2, 5., 0.05)
x[1] = np.pi/2. # eliptica

Bm = []
for x[0] in Rs:
    B = bp.return_B(x)
    Bm += [ np.square(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]) ]

Bm = np.array(Bm)
"""
#EOF
