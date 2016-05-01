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


#--- set B-turbulence model
pd['n_modos']    = 128
pd['lambda_min'] = ((5e-5)*AU_in_cm)
#--- corregimos input
psim['atol']     = 1e-5 #3.3e-6 #1e-6 #1e-5
psim['rtol']     = 0.0 #1e-6
psim['tmax']     = 1e4 #0.3e4 #4e4
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

    print " [r:%d] now i'll write" % rank

    fo['mine'] = 1    # file ocupado
    dpath = 'pla%03d/' % npla
    ff.SaveToFile(m, dpath=dpath, f=fo)
    pause(0.1)
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
