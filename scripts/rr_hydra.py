#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import src.cython_wrapper as cw
import shared.funcs as sf
import mhandlers as mh
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

#++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++
# IMPORTANT NOTES:
# * If you change the dataset structure of the output
#   file (ie. the HDF5 tree), then remember to also 
#   change the 'version' field inside! 
#++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++


#--- globals
AUincm = 1.5e13        # [cm]
r2d    = 180./np.pi    # rad to degrees


#--- parse args
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
    Trace particle trajectories in different magnetic 
    field scenarios/configurations.
    """,
)
# subparsers
subparsers = parser.add_subparsers(
description="""
Use one of the submodules below.
""",
)
# config the subparsers
for mod_name in ['slab2d', ]:
    # grab the class
    mod = getattr(mh, 'dBmodel__'+mod_name)()
    # grab the parser of that class
    subparser_ = subparsers.add_parser(
    mod_name, 
    help=mod.help,
    parents=[mod.parser], 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # all subparser will have this options in common:
    subparser_.add_argument(
    '-Nth',
    type=int,
    default=0, #16,
    help='number of theta values for the initial velocities.',
    )
    subparser_.add_argument(
    '-tho',
    type=float,
    default=None,
    help='theta angle for a single pitch angle. If used, the -Nth cannot be used.',
    )
    subparser_.add_argument(
    '-Nph',
    type=int,
    default=8,
    help='number of phi values for the initial velocities.',
    )
    subparser_.add_argument(
    '-tmax',
    type=float,
    default=4e4,
    help='max time for simulation',
    )
    subparser_.add_argument(
    '-tscale',
    default='mixed', # 'linear',
    help='time scale for the saved trajectories; can be mixed or linear.',
    )
    subparser_.add_argument(
    '-nsave',
    type=int,
    default=20,
    help="""For -tscale=linear, the total number of saved times; for 
    -tscale=mixed, the number of saved times in each decade.""",
    )
    subparser_.add_argument(
    '-eps',
    type=float,
    default=4.64e-6, # 3.33e-6
    help='epsilon parameter; i.e. (error-step)/(lambda_min)',
    )
    subparser_.add_argument(
    '-trails',
    action='store_true',
    default=False,
    help='whether or not to grab the trails data; that is, trajectory \
    traces of the particles just before they made a back-scatter (i.e. \
    when the particle changes sign in pitch angle)',
    )
    subparser_.add_argument(
    '-taubd',
    type=float,
    nargs='+',
    default=None, #[0.8, 1.1, 1.4, 1.6],
    help='list of borders for the bands of interest in the collision-time histogram. \
    This works if the code was compiled with the WATCH_TRAIL preprocessor. \
    If used, the code will collect particle trails (positions along its trajectory) when \
    it backscatters and the associated collision-time is within these bands of interest. \
    The values of this bands are assumed to be in units of gyro-cycles.'
    )
    subparser_.add_argument(
    '-of', '--fname_out',
    type=str,
    default='../out/test/testt.h5',
    help='output HDF5 file.',
    )

pa = parser.parse_args()

# grab the name of the module selected in the CLI
mod_selected = sys.argv[1]
print "\n ---> Using %s\n"%mod_selected 
mod = getattr(mh, 'dBmodel__'+mod_selected)()

#--- orientations of the initial velocities
th, ph = sf.generate_init_conds(tho=pa.tho/r2d, Nth=pa.Nth, Nph=pa.Nph)
mu = np.cos(th)         # pitch-angle cosine

#--- MPI
#comm        = MPI.COMM_WORLD
#rank        = comm.Get_rank()   # proc rank
#wsize       = comm.Get_size()   # number of proc
#rank  = int(sys.argv[1])
#wsize = int(sys.argv[2])

# build the times series
getattr(mod, 'run')(pa, pd=pd, psim=psim, 
    pother=pother, 
    mu=mu, 
    ph=ph,
    comm=MPI.COMM_WORLD,
    nB=nB)



#EOF
