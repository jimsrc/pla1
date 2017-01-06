"""
parameters file
"""
import numpy as np
from os import environ

AU_in_cm = 1.5e13

fname_ori = '{PLA1}/scripts/orientations_122.in'.format(**environ)
#ori = np.loadtxt('orientations_isotropic_Nth16_Nph8.in')
ori = np.loadtxt(fname_ori)
mu  = np.cos(ori[:,0])
ph  = ori[:,1]

#--- B-field parameters only
pd = {
    # parametros fisicos
    'Nm_slab'       : 128,
    'Nm_2d'         : 128,
    'lmin_s'        : ((5e-5)*AU_in_cm),
    'lmax_s'        : 1.0*AU_in_cm,
    'lmin_2d'       : ((5e-5)*AU_in_cm),
    'lmax_2d'       : 1.0*AU_in_cm,
    'Lc_slab'       : 0.01*AU_in_cm,
    #'Lc_2d'         : 0.01*AU_in_cm,
    'sigma_Bo_ratio': 1.0,
    #'percent_slab'  : 0.2,
    #'percent_2d'    : 0.8,
    #'Bo'            : 5e-5,   # [Gauss]
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
    #'rigidity'      : 4.44583E+08,  # [V]
    'tmax'          : 4e4,
    'FracGyroperiod': 5e-2,
    'hmin'          : 0.0,
    'mu'            : mu[npla],
    #'ph'            : ph[npla],
    'rtol'          : 1e-5,
    'atol'          : 1e-5,
}

pother = {
    'str_timescale' : 'mixed',
    'nsave'         : 20,
    'tmaxHistTau'   : 150,
    'nHist'         : 150,
    'nThColl'       : 100,
    'i'             : npla,
    'j'             : nB,
    'dir_out'       : './',
}
