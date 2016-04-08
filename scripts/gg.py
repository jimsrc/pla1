#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from h5py import File as f5
import os
from os.path import isfile, isdir

class dummy:
    pass

def read_pla(f=None, iB=0, ipla=0):
    if f==None:
        print " no file dude!"
        raise SystemExit

    path = 'B%02d/pla%03d' % (iB, ipla)
    if not(hasattr(f, path)):
        print " ---> file doesn't have '%s' attribute!!" % path
        raise SystemExit

    v     = dummy()
    v.t   = f[path+'/t'][...]
    v.xyz = f[path+'/xyz'][...].T
    v.mu  = f[path+'/mu'][...]
    v.err = f[path+'/err'][...]
    v.ntot_B   = f['ntot_B'].value
    v.ntot_pla = f['ntot_pla'].value
    #return t, mu, err
    return v
"""  -----> revisar con:
/home/jim/simulacion/pla_stochastic/composite_turbulence/massive_orbits/output/Ek.1.0e+08eV/Nm128/slab0.20/sig.1.0e+00/Lc2D.1.0e-02_LcSlab.1.0e-02/B00_pla083_traj.dat

"""


PLAS     = os.environ['PLAS']
Ek       = 1e8 #1e8 #1e7 #1e6
slab     = 0.2
Nm       = 128
sig      = 1.0
Lc2d, LcSlab = 1e-2, 1e-2

name_Ek = 'Ek.%1.1eeV' % Ek
name_Nm = 'Nm%03d' % Nm
name_slab = 'slab%2.2f' % slab
name_sig  = 'sig.%1.1e' % sig
name_Lc   = 'Lc2D.%1.1e_LcSlab.%1.1e' % (Lc2d, LcSlab)

dir_src  = '%s/output' % PLAS + '/' + name_Ek
dir_src  += '/' + name_Nm + '/' + name_slab + '/' + name_sig + '/' + name_Lc; d2 = dir_src
fname_inp_base = '%s__%s__%s__%s__%s.h5' % (name_Ek, name_Nm, name_slab, name_sig, name_Lc)
fname_inp = '%s/%s' % (dir_src, fname_inp_base)
#fname_inp = '%s/out.h5' % (dir_src)

#
dir_src = '%s/output/Ek.1.0e+06eV_rtol.1e-05' % PLAS
fname_inp = '%s/out.h5' % dir_src
#

f = f5(fname_inp, 'r')
#t, mu, err = read_pla(f, iB=0, ipla=0)
#v = read_pla(f, iB=0, ipla=0)

#EOF
