#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
import funcs as ff
from pylab import pause
#--- parameters
from params import (
    nB, pd, psim, pother, mu, ph, AU_in_cm
)
from mpi4py import MPI
from h5py import File as h5
from os.path import isfile, isdir
import os


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
psim['rigidity'] = 4.33306E+07
psim['tmax']     = 4e4 #0.3e4 #4e4
rl = cw.calc_Rlarmor(psim['rigidity'],pd['Bo']) #[cm]
eps_o = 3.3e-6 #3.33e-5 #4e-5 # ratio: (error-step)/(lambda_min)
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
comm        = MPI.COMM_WORLD
rank        = comm.Get_rank()   # proc rank
wsize       = comm.Get_size()   # number of proc
#---------------------------
if rank==0:
    print " simul parameters\n", psim

pla_bd  = ff.equi_bounds(0, mu.size-1, wsize) # bounds
plas    = np.arange(pla_bd[rank], pla_bd[rank+1]) # plas for each proc

if rank==0 and isfile(fname_out):  # backup if already exists
    os.system('mv {fname} {fname}_'.format(fname=fname_out))
#new = True # new-file flag

for npla in plas: #[25:]:
    #--- set particle id && direction
    pother['i']  = npla
    psim['mu']   = mu[npla]
    psim['phi']  = ph[npla]

    m.build(**pother)

    m.SetSim(**psim)
    m.RunSim()

    print " [r:%d] simulation (pla:%d) finished!" % (rank, npla)

    if isfile(fname_out):
        pause(0.2*(rank + 1))
        fo = h5(fname_out, 'r+') # append
    else:
        fo = h5(fname_out, 'w')  # create

    while 'mine' in fo.keys():
        pause(0.1)
        fo = h5(fname_out, 'r+') # read again

    fo['mine'] = 1    # file ocupado
    print " [r:%d] now i'll write" % rank
    dpath = 'pla%03d/' % npla
    ff.SaveToFile(m, dpath=dpath, f=fo)
    #--- histos of step-size
    st      = m.step_save # shape: (2,Odeint::MAXSTP)
    st_tot  = st[0,:][st[0,:]!=0] # filter-out zero values
    st_part = st[1,:][st[1,:]!=0] # ""
    #st     = st[0,:][st!=0]   # filter-out zero values
    nbin = 1000
    h    = np.histogram2d(st_tot, st_part, 
                bins=[nbin, nbin], 
                normed=False)
    HStp = h[0].T # filter-out zero values
    bins_StepTot  = 0.5*(h[1][:-1] + h[1][1:])
    bins_StepPart = 0.5*(h[2][:-1] + h[2][1:])
    fo[dpath+'HistStep/HStep'] = HStp.sum(axis=0)
    #fo[dpath+'HistStep/bins_StepTot']  = bins_StepTot
    fo[dpath+'HistStep/bins_StepPart'] = bins_StepPart
    fo[dpath+'HistStep/nbin']          = nbin
    pause(0.3)
    fo.pop('mine')    # libera file

    fo.close()
    print " [r:{rank}] closing: {fname}".format(rank=rank,fname=fname_out)

    # everybody waits for each processor to finish the i/o
    #comm.Barrier() 

print " [r:%d] I'm finished!" % rank
# clean backup
if rank==0:
    os.system('rm {fname}_'.format(fname=fname_out))
#EOF
