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

"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01
1e10    1.0898E+10   4.853544E-02   9.963142E-01
"""
#--- set B-turbulence model
pd['n_modos']    = 128
pd['lambda_min'] = ((5e-5)*AU_in_cm)
#--- corregimos input
psim['rigidity'] = 1.69604E+09
psim['tmax']     = 4e4 #0.3e4 #4e4
rl = cw.calc_Rlarmor(psim['rigidity'],pd['Bo']) #[cm]
eps_o = 1.0e-5 #3.33e-6 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
psim['atol']     = pd['lambda_min']*eps_o/rl
psim['rtol']     = 0.0 #1e-6
#---------------------------
#--- output
po = {}
po.update(psim)
po.update(pd)
po['lambda_min'] /= AU_in_cm
dir_out = '.'
fname_out = dir_out+'/R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_Nm.{n_modos:04d}_lmin.{lambda_min:1.1e}.h5'.format(**po)

# call simulator
m = cw.mgr()
m.set_Bmodel(pdict=pd, nB=nB)

#--- MPI
#comm        = MPI.COMM_WORLD
#rank        = comm.Get_rank()   # proc rank
#wsize       = comm.Get_size()   # number of proc
rank = int(sys.argv[1])
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
