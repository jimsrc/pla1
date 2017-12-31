#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from h5py import File as h5
from os.path import isfile, isdir
import numpy as np
from numpy import power as pow
from pylab import (
    pause, find, figure, close
)
from datetime import datetime, timedelta
#from numpy import min, max
import os, sys
from glob import glob
from Bparker.Bparker import return_B as Bparker_vector
from numpy.linalg import norm


def memory_stat(SNUMB, dtype=np.float64, wsize=1):
    """
    routine to check whether we have enough available RAM for
    the data we are asking to allocate.
    NOTE that we are estimating these in terms of those
    pieces of C++ code where there is a "WATCH_MEMORY" comment (see
    .cc sources!)
    """
    #--- memory check
    # NOTE:
    # 1024 bytes = 1 KB
    # 1024 KB    = 1 MB
    # 1024 MB    = 1 GB
    # 1GB        = 1024*1024*1024 bytes
    SDOUB    = np.dtype(dtype).itemsize  # [bytes] size of a double (float64 in python)
    TOTBYTES = SNUMB*SDOUB              # [bytes] total size per processor
    # NOTE: we are multiplying by the **total number of nodes** 'wsize'
    STOT     = wsize*TOTBYTES/(1024.*1024.*1024.)  # [GB] total size of our data
    #--- check with available RAM
    import psutil
    avail_ram = psutil.virtual_memory()[1]/(1024.*1024.*1024.)

    # True: if we have more available RAM than the
    #       requested space.
    return STOT, avail_ram


def calc_Rlarmor(rigidity, Bo, full=False):
    """
    input:
    Ek      : [eV] kinetic energy
    rigi..  : [V] rigidity
    Bo      : [G] magnetic field in Gauss
    output:
    Rl  : [cm] larmor radii
    """
    q = (4.8032*1e-10) # [statC] carga PROTON
    mo = 1.6726e-24 # [gr] masa PROTON
    c = 3e10            # [cm/s] light speed
    AU_in_cm = 1.5e13     # [cm]
    E_reposo=938272013.0  # [eV] PROTON
    #beta, gamma, omg, v

    #rigidity = sqrt(Ek*Ek + 2.*Ek*E_reposo);
    #------------------------CALCULO DE GAMMA Y BETA
    gamma = pow(pow(rigidity/E_reposo,2) + 1. , 0.5)
    beta = pow(1. - 1/(gamma*gamma) , 0.5)
    #------------------------------CALCULO CICLOTRON
    omg = q * Bo / (gamma * mo * c)     # [s^-1]
    #---------------------------CALCULO RADIO LARMOR
    v   = beta * c              # [cm/s]
    #Rl[0]  = (v / omg) /AU_in_cm  # [AU]
    if full:
        return (v/omg), v, omg
    else:
        return (v / omg) # [cm]

def Lc_memilia(r=1.0):
    """ 
    formula from Maria Emilia's thesis (Sec 4.4, p.88).
    """
    Lc = 0.89*(r**(0.43))*1e6/(1e8) # [AU]
    return Lc


def Bo_parker(r=1.0, th=np.pi/2., ph=0.0):
    """ input:
    r  [AU]  : heliodistance
    th [rad] : spherical polar angle (co-latitude)
    ph [rad] : azimuth
    default: (r=1.0, th=pi/2, ph=0.0)
    """
    #--- calc B-field at position (r, th)
    pos = np.array([r, th, ph], dtype=np.float32) #[AU], [rad], [rad]: position
    Bo  = norm(Bparker_vector(pos)) # [Gauss]
    return Bo


def unify_all(fname_out, psim, wsize):
    """ function to unify all the
        output .h5 of all processors.
	NOTE: this is for working in
	HYDRA cluster only (it has no -fPIC MPI-configured).
    """
    fout = h5(fname_out, 'w')
    for r in range(wsize):
        fnm_inp = fname_out+'_%02d'%r
        finp = h5(fnm_inp, 'r')
        cont = finp.keys()       # list of groups
        for c in cont:           # iterate over each group
            finp.copy(c, fout)
        finp.close()
        print(" [*] removing partial file: "+fnm_inp)
        os.system('rm {fname}'.format(fname=fnm_inp))

    for pnm in psim.keys():
        fout['psim/%s'%pnm] = psim[pnm]

    print "\n [+] We generated: "+fout.filename
    fout.close()
    _right_now = datetime.now().strftime("%d %b %Y %H:%M:%S")
    print "\n ----> FINISHED UNIFYING OUTPUT @ " + _right_now +'\n'


def equi_bounds(dini, dend, n):
    """
    bounds in 'n' equipartitioned segments inside a numerical
    range between 'dini' and 'dend'; starting in 'start'
    """
    interv = equi_intervals(dini, dend, n).cumsum()
    bd = []
    bd += [ dini ]
    for it in interv:
        bd += [ it ] 

    bd = np.array(bd)
    return bd


def equi_intervals(dini, dend, n): 
    """ 
    returns most equi-partitioned tuple in the numerical
    range between 'dini' and 'dend'
    """
    #days = (dend - dini).days
    days = dend - dini + 1
    days_part = np.zeros(n, dtype=np.int)
    resid = np.mod(days, n)
    for i in range(n-resid):
        days_part[i] = days/n

    # last positions where I put residuals
    last = np.arange(start=-1,stop=-resid-1,step=-1)
    for i in last:
        days_part[i] = days/n+1

    assert np.sum(days_part)==days, \
        " --> somethng went wrong!  :/ "

    return days_part


def w2file(rank, fname, var, vname, mode='r+'):
    try:
        pause(0.2*(rank+1))
        f = h5(fname, mode=mode)
        f[vname] = var
        f.close()
        #print " [r:%d] ########## TERMINE ########## " % rank
    except NameError:
        pause((rank+1)*0.3)
        f.close()
        w2file(rank, fname, var, vname, mode='r+')
    except RuntimeError as e:
        if e.message.startswith("Unable to create link"):
            pass
        pause((rank+1)*0.3)
        if isfile(fname):
            w2file(rank, fname, var, vname, mode='r+')
        else:
            w2file(rank, fname, var, vname, mode='w')
    except IOError as e:
        #print " [r:%d] MSG: %s" % (rank, e.message)
        pause((rank+1)*0.3)
        #print " [r:%d] --->"%rank, dir(e)
        if isfile(fname):
            w2file(rank, fname, var, vname, mode='r+')
        else:
            w2file(rank, fname, var, vname, mode='w')


def SaveToFile_ii(rank, m, dpath, fname_out):
    vsave = {
        'tadim'     : m.tsave,      # [omega^-1]
        'ysave'     : m.ysave,      # [1]
        'xyz'       : m.xyz,        # [1]
        'err'       : m.err,        # [1] no porcentaje
        'mu'        : m.mu,         # [1]
        'vel'       : m.vel,        # [scl_vel] normalizadas
        'scl_wc'    : m.scl['wc'],  # [s^-1]
        'scl_vel'   : m.scl['vel'], # [cm/s]
        'scl_beta'  : m.scl['beta'], # [1]
        'scl_gamma' : m.scl['gamma'], # [1]
        'scl_rl'    : m.scl['rl'],  # [cm]
    }

    #print " ---> saving: " + f.filename
    for name in vsave.keys():
        dname = dpath + name   # anteponemos algo
        w2file(rank, fname_out, vsave[name], dname)
        #f[dname] = vsave[name]


def SaveToFile(m, dpath='', f=None, nbin_step=None, btrails=False):
    """
    input:
    - m         : cython-wrapper-module
    - dpath     : directory-path for .h5 file
    - f         : output file object
    - nbin_step : number of bins for step-size histos
    - btrails   : whether or not to save trails data
    """
    vsave = {
    'tadim'     : m.tsave,      # [omega^-1]
    'xyz'       : m.xyz,        # [1]
    'err'       : m.err,        # [1] no porcentaje
    'mu'        : m.mu,         # [1]
    'vel'       : m.vel,        # [scl_vel] normalizadas
    'trun'      : m.trun,       # [sec] runtime
    }

    #--- variables
    for name in vsave.keys():
        dname = dpath + name   # anteponemos algo
        f[dname] = vsave[name]

    #--- histos of step-size
    st      = m.step_save # shape: (2,Odeint::MAXSTP)
    st_tot  = st[0,:][st[0,:]!=0] # filter-out zero values
    st_part = st[1,:][st[1,:]!=0] # ""
    h       = np.histogram2d(st_tot, st_part, 
              bins=[nbin_step, nbin_step], 
              normed=False)
    HStp          = h[0].T # filter-out zero values
    bins_StepTot  = 0.5*(h[1][:-1] + h[1][1:])
    bins_StepPart = 0.5*(h[2][:-1] + h[2][1:])
    f[dpath+'HistStep/HStep'] = HStp.sum(axis=0)
    f[dpath+'HistStep/bins_StepPart'] = bins_StepPart
    f[dpath+'HistStep/nbin'] = nbin_step

    #--- histos for tau-collision
    tauLg   = np.log10(m.Tau[:,1]) # log([1/omega?])
    h       = np.histogram(tauLg, bins=nbin_step, normed=False)
    hc      = h[0]
    hbin    = 0.5*(h[1][:-1] + h[1][1:])    # log([1/omega?])
    f[dpath+'HistTau_log'] = np.array([hc, hbin])

    #--- histos theta (angle between x-y plane and z-axis, @collision)
    theta   = m.Tau[:,4] # [deg]
    h       = np.histogram(theta, bins=nbin_step, normed=False)
    hc      = h[0]
    hbin    = 0.5*(h[1][:-1] + h[1][1:])    # [deg]
    f[dpath+'HistThetaColl'] = np.array([hc, hbin])

    #--- grab the trails if the code had it activated
    if btrails:
        print m.ptrails.shape
        import pdb; pdb.set_trace()



def save_to_h5(m, fname=None, dpath='', file=None, close=True):
    vsave = {
        'tadim'     : m.tsave,      # [omega^-1]
        'ysave'     : m.ysave,      # [1]
        'xyz'       : m.xyz,        # [1]
        'err'       : m.err,        # [1] no porcentaje
        'mu'        : m.mu,         # [1]
        'vel'       : m.vel,        # [scl_vel] normalizadas
        'scl_wc'    : m.scl['wc'],  # [s^-1]
        'scl_vel'   : m.scl['vel'], # [cm/s]
        'scl_beta'  : m.scl['beta'], # [1]
        'scl_gamma' : m.scl['gamma'], # [1]
        'scl_rl'    : m.scl['rl'],  # [cm]
    }
    if isfile(fname):
        print " ---> Overwriting: " + fname 

    if file is None:
        f = h5(fname, 'w')
    else:
        f = file

    print " ---> saving: " + fname 
    for name in vsave.keys():
        name = dpath + name   # anteponemos algo
        f[name] = vsave[name]

    if close:
        f.close()
        print " ---> file closed: ", fname
        return 0
    else:
        return f

#EOF
