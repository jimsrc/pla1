#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
import shared.funcs as ff
from pylab import pause
#--- parameters
from params import (
    nB, pd, psim, pother, AU_in_cm
)
from mpi4py import MPI
from h5py import File as h5
from os.path import isfile, isdir
import os, sys, argparse
from glob import glob
from numpy import power, log10

#--- globals
AUincm = 1.5e13             # [cm]


#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-of', '--fname_out',
type=str,
default='../out/test/testt.h5',
help='output HDF5 file.',
)
parser.add_argument(
'-NmS', '--Nm_slab',
type=int,
default=128,
help='number of slab modes',
)
parser.add_argument(
'-Nm2d', '--Nm_2d',
type=int,
default=128,
help='number of 2D modes',
)
parser.add_argument(
'-Nth', '--Nth',
type=int,
default=16,
help='number of theta values for the initial velocities.',
)
parser.add_argument(
'-Nph', '--Nph',
type=int,
default=8,
help='number of phi values for the initial velocities.',
)
parser.add_argument(
'-ro', '--ro',
type=float,
default=0.9,
help='heliodistance in AU.',
)
parser.add_argument(
'-tmax', '--tmax',
type=float,
default=4e4,
help='max time for simulation',
)
parser.add_argument(
'-eps', '--eps',
type=float,
default=4.64e-6, # 3.33e-6
help='epsilon parameter; i.e. (error-step)/(lambda_min)',
)
parser.add_argument(
'-trails', '--trails',
action='store_true',
default=False,
help='whether or not to grab the trails data; that is, trajectory \
traces of the particles just before they made a back-scatter (i.e. \
when the particle changes sign in pitch angle)',
)
pa = parser.parse_args()


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

lc = ff.Lc_memilia(r=pa.ro)   # [AU], this gives the correlation-LENGTH
"""
from result in '16d5869' commit (from PLAS repo), we have the relation:
    y = m*x + b
    where,
    x = log10(lambda_c)  # REAL/FITTED correlation length
    y = log10(Lc)        # correlation scale
Then:
    Lc = 10.**(m*log10(lambda_c) + b)
    with:
    m = 1.026
    b = 0.2755
NOTE: this is valid in the range log10(Lc):[-2.5, 1.0]
"""
# obtain the correlation-SCALE as function of the correlation-LENGTH
Lc_slab = power(10., 1.026*log10(lc) + 0.2755) # from the above comments

Rl = cw.calc_Rlarmor(
    rigidity=1.69604E+09, #1.69604E+09,    # [V]
    Bo=ff.Bo_parker(r=pa.ro)    # [Gauss]
    )/AUincm                 # [AU] Larmor radii
#--- set B-turbulence model
pd.update({
'Nm_slab'       : 128,
'Nm_2d'         : 128,
'lmin_s'        : 5e-5/Rl, #[lmin_s/Rl] 
'lmax_s'        : pa.ro/Rl,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-5/Rl, #[lmin_2d/Rl] 
'lmax_2d'       : pa.ro/Rl,  #[lmax_2d/Rl] 
'Lc_slab'       : Lc_slab/Rl,  # in units of Larmor-radii
'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
'sigma_Bo_ratio': 0.3, # [1] fluctuation energy
'ratio_slab'    : 0.2, # [1] (energy_slab)/(energy_total)
})
#--- corregimos input
psim['tmax']     = pa.tmax     # [1/omega]
eps_o            = pa.eps      # ratio: (error-step)/(lambda_min)
lmin             = np.min([pd['lmin_s'], pd['lmin_2d']]) # [cm] smallest turb scale
psim['atol']     = lmin*eps_o  # [1]
psim['rtol']     = 0.0 #1e-6


#--- orientations of the initial velocities
ph, th = ff.generate_init_conds(Nth=pa.Nth, Nph=pa.Nph)
mu = np.cos(th)         # pitch-angle


#--- output
po = {}
po.update(psim)
po.update(pd)
# add some stuff
po.update({
'r'     : pa.ro,           # [AU] heliodistance
'eps_o' : eps_o,        # [1]  precision
'lmin'  : lmin,         # [1]  minimum turb scale (*)
'RloLc' : Rl/Lc_slab,   # [1]  (r_larmor)/(Lc_slab)
})
# (*) according to the definitions above of 'lmin_*' they are
# adimensional but they're in units of Rl.

odir = '/'.join(pa.fname_out.split('/')[:-1])
assert isdir(odir), ' [-] doesnt exist output directory:\n %s\n'%odir

#--- call simulator
m = cw.mgr()
m.set_Bmodel(pdict=pd, nB=nB)

#--- MPI
comm        = MPI.COMM_WORLD
rank        = comm.Get_rank()   # proc rank
wsize       = comm.Get_size()   # number of proc
#rank  = int(sys.argv[1])
#wsize = int(sys.argv[2])


#--- memory check
# NOTE: we estimate the "requested data" in terms
# of those pieces of C++ code where there is a
# "WATCH_MEMORY" comment (see .cc sources!).
STOT, avail_ram = ff.memory_stat(
    SNUMB=(2)*m.MAXSTP, # according pieces in "CHECK_MEMORY" comments
    dtype=np.float64,   # np.float64==double 
    wsize=wsize)
if rank==0:
    print """
    +++++++++++ memory report +++++++++++++
     You are asking for : {mem_ask} GB (estimated of data)
     but we have        : {mem_avail} GB (available RAM)
    +++++++++++++++++++++++++++++++++++++++
    """.format(mem_ask=STOT, mem_avail=avail_ram)
    if (STOT>avail_ram):
        print "\n --> ERROR: NOT ENOUGH MEMORY!\n"
        # not the most correct to kill it, but works
        comm.Abort()
# wait for the above check
comm.Barrier()


#---------------------------
if rank==0:
    print " simul parameters\n", psim

pla_bd  = ff.equi_bounds(0, mu.size-1, wsize) # bounds
plas    = np.arange(pla_bd[rank], pla_bd[rank+1]) # plas for each proc

if rank==0 and isfile(pa.fname_out):  # backup if already exists
    os.system('mv {fname} {fname}_'.format(fname=pa.fname_out))

fname_out_tmp = pa.fname_out+'_%02d'%rank # output of each processor
fo = h5(fname_out_tmp, 'w') 
nbin = 1000 # for step-size histograms
for npla in plas: #[25:]:
    #--- set particle id && direction
    pother['i']  = npla
    psim['mu']   = mu[npla]
    psim['ph']   = ph[npla]

    m.build(**pother)

    m.SetSim(**psim)
    m.RunSim()

    print " [r:%03d] (pla:%d of %d) finished!" % (rank, npla, plas[-1])
    dpath = 'pla%03d/' % npla
    print " [r:%03d] writing: %s" % (rank, fo.filename)
    ff.SaveToFile(m, dpath, fo, nbin_step=nbin, btrails=pa.trails)
    m.clean()

print " [r:%03d] closing: %s" % (rank, fo.filename)
fo.close()
pause(1)
os.system('touch %s_finished'%fname_out_tmp)
print " [r:%03d] I'm finished!" % rank

if rank==0:
    #--- let's check if everyone has already finished
    n = 0
    while n<wsize:
        pause(2) # wait before next check
        fs = glob(pa.fname_out+'_*_finished') #list of files
        n = len(fs)
    
    #--- now proceed to unify
    ff.unify_all(pa.fname_out, psim=po, wsize=wsize)
    #--- now we can delete the "*_finished" files
    for fnm in fs:
        os.system('rm '+fnm)

#EOF
