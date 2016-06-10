#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
import shared.funcs as ff
from pylab import pause
#--- parameters
from params import (
    nB, pd, psim, pother, mu, ph, AU_in_cm
)
#from mpi4py import MPI
from h5py import File as h5
from os.path import isfile, isdir
import os, sys
from glob import glob
#from Bparker.Bparker import return_B as Bparker

"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01  (*)
1e10    1.0898E+10   4.853544E-02   9.963142E-01
--- Bparker model
(Ek=1e9eV)
r[AU]    B[nT]       Rl[AU]         Lc[AU]      Rl/Lc   Rl/(5e-5AU)
0.2      91.0542857  4.147812E-04   0.0044549   0.09    8.30
0.3      41.4298571  9.116035E-04   0.0053034   0.17    18.23
0.4      24.041      1.570966E-03   0.0060017   0.26    31.42
0.5      15.972      2.364613E-03   0.0066061   0.36    47.29
0.7      8.897       4.244982E-03   0.0076345   0.56    84.90
0.8      7.14642857  5.284822E-03   0.0080857   0.65    105.70
1.0      5.0         7.553521E-03   0.0089      0.85    151.07
2.0      1.99653571  1.891657E-02   0.0119904   1.58    378.33
"""
#--- set B-turbulence model
pd.update({
    'Nm_slab'       : 128,
    'Nm_2d'         : 128,
    'lmin_s'        : (5e-5)*AU_in_cm,
    'lmax_s'        : (1.0) *AU_in_cm,
    'lmin_2d'       : (5e-5)*AU_in_cm,
    'lmax_2d'       : (1.0) *AU_in_cm,
    'Lc_slab'       : 0.01*AU_in_cm,
    'Lc_2d'         : 0.01*AU_in_cm,
    'sigma_Bo_ratio': 1.0,
    'Bo'            : ff.Bo_parker(r=0.3),   # [Gauss]
})
#--- corregimos input
psim['rigidity'] = 1.69604E+09
psim['tmax']     = 4e4 #0.3e4 #4e4
rl = cw.calc_Rlarmor(psim['rigidity'],pd['Bo']) #[cm]
eps_o = 1.0e-5 #3.33e-6 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
lmin             = np.min([pd['lmin_s'], pd['lmin_2d']]) # smallest turb scale
psim['atol']     = lmin*eps_o/rl
psim['rtol']     = 0.0 #1e-6
#---------------------------
#--- output
po = {}
po.update(psim)
po.update(pd)
po['lmin_s'] /= AU_in_cm
dir_out = '.'
fname_out = dir_out+'/R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_NmS.{Nm_slab:04d}_Nm2d.{Nm_2d:04d}_lminS.{lmin_s:1.1e}.h5'.format(**po)

# call simulator
m = cw.mgr()
m.set_Bmodel(pdict=pd, nB=nB)

#--- MPI
#comm        = MPI.COMM_WORLD
#rank        = comm.Get_rank()   # proc rank
#wsize       = comm.Get_size()   # number of proc
rank  = int(sys.argv[1])
wsize = int(sys.argv[2])
#---------------------------
if rank==0:
    print " simul parameters\n", psim

pla_bd  = ff.equi_bounds(0, mu.size-1, wsize) # bounds
plas    = np.arange(pla_bd[rank], pla_bd[rank+1]) # plas for each proc

if rank==0 and isfile(fname_out):  # backup if already exists
    os.system('mv {fname} {fname}_'.format(fname=fname_out))

fname_out_tmp = fname_out+'_%02d'%rank # output of each processor
fo = h5(fname_out_tmp, 'w') 
nbin = 1000 # for step-size histograms
for npla in plas: #[25:]:
    #--- set particle id && direction
    pother['i']  = npla
    psim['mu']   = mu[npla]
    psim['phi']  = ph[npla]

    m.build(**pother)

    m.SetSim(**psim)
    m.RunSim()

    print " [r:%d] simulation (pla:%d) finished!" % (rank, npla)
    dpath = 'pla%03d/' % npla
    print " [r:%d] now i'll write" % rank
    ff.SaveToFile(m, dpath, fo)

    #--- histos of step-size
    st      = m.step_save # shape: (2,Odeint::MAXSTP)
    st_tot  = st[0,:][st[0,:]!=0] # filter-out zero values
    st_part = st[1,:][st[1,:]!=0] # ""
    h       = np.histogram2d(st_tot, st_part, 
              bins=[nbin, nbin], 
              normed=False)
    HStp          = h[0].T # filter-out zero values
    bins_StepTot  = 0.5*(h[1][:-1] + h[1][1:])
    bins_StepPart = 0.5*(h[2][:-1] + h[2][1:])
    fo[dpath+'HistStep/HStep'] = HStp.sum(axis=0)
    fo[dpath+'HistStep/bins_StepPart'] = bins_StepPart
    fo[dpath+'HistStep/nbin'] = nbin
    print " [r:{rank}] closing: {fname}".format(rank=rank,fname=fname_out)

fo.close()
pause(1)
os.system('touch %s_finished'%fname_out_tmp)
print " [r:%d] I'm finished!" % rank

if rank==0:
    #--- let's check if everyone has already finished
    n = 0
    while n<wsize:
        pause(2) # wait before next check
        fs = glob(fname_out+'_*_finished') #list of files
        n = len(fs)
    
    #--- now proceed to unify
    ff.unify_all(fname_out, wsize)
    #--- now we can delete the "*_finished" files
    for fnm in fs:
        os.system('rm '+fnm)

#EOF
