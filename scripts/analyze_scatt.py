#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import figure, close
import h5py, argparse


# retrieve the IDentifiers :-)
parser = argparse.ArgumentParser(
description="""
analyze the trails and scatterings associated.
""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-fi', '--fname_inp',  # the .h5 file is assumed to have a '%03d' format for the ID
type=str,
help='input HDF5 file.',
)
pa = parser.parse_args()



pi2 = 2.0*np.pi
f = h5py.File(pa.fname_inp, 'r')

# number of pla, and number of bounce
ipla, ireb  = 0, 5

tr          = f['pla%03d/trails/trajecs'%ipla][...]
tau_b       = f['pla%03d/trails/tau_b'%ipla][...]
t, mu, xyz  = tr[ireb,:,0], tr[ireb,:,1], tr[ireb,:,2:]
x, y, z     = xyz.T
#cc          = tau_b > (0.8*pi2)
ccmu        = (t[0]-t) < tau_b[ireb]
tt          = (t[0] - t)/pi2            # units of gyrocycles respect to the last collision

mu_avr, mu_std = np.mean(mu), np.std(mu)

fig     = figure(1, figsize=(6,4))
ax      = fig.add_subplot(111)

ax.plot(tt, mu, '-')
ax.plot(tt[ccmu], mu[ccmu], 'r-')
ax.axhline(y=0.0, lw=3., c='k', alpha=0.6)

ax.grid()

fig.savefig('./test.png', bbox_inches='tight')
close(fig)

#EOF
