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
pd.update({
    'Nm_slab'       : 128,
    'Nm_2d'         : 128,
    'lmin_s'        : ((5e-5)*AU_in_cm),
    'lmax_s'        : 1.0*AU_in_cm,
    'lmin_2d'       : ((5e-5)*AU_in_cm),
    'lmax_2d'       : 1.0*AU_in_cm,
})

#--- corregimos input
psim['rigidity'] = 4.33306E+07
psim['tmax']     = 4.4e4 #0.3e4 #4e4
rl = cw.calc_Rlarmor(psim['rigidity'],pd['Bo']) #[cm]
eps_o = 3.33e-4 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
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
fname_out = dir_out+'/R.{rigidity:1.2e}_atol.{atol:1.1e}_rtol.{rtol:1.1e}_NmS.{Nm_slab:04d}_lminS.{lmin_s:1.1e}.h5'.format(**po)

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
    #ff.SaveToFile_ii(rank, m, dpath, fname_out)

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
print " [r:%d] I'm finished!" % rank
comm.Barrier() # I need everyone to stop touching the files
# let's unify all files into one
if rank==0:
    fout = h5(fname_out, 'w')
    for r in range(wsize):
        fnm_inp = fname_out+'_%02d'%r
        finp = h5(fnm_inp, 'r')
        cont = finp.keys()       # list of groups
        for c in cont:           # iterate over each group
            finp.copy(c, fout)
        finp.close()
        os.system('rm {fname}'.format(fname=fnm_inp))

    print " ----> We generated: "+fout.filename
    fout.close()
    # clean backup
    os.system('rm {fname}_'.format(fname=fname_out))
    print " [r:%d] I'm finished!" % rank
#EOF
