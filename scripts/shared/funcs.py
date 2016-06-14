#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from h5py import File as h5
from os.path import isfile, isdir
import numpy as np
from pylab import pause, find
#from numpy import min, max
import os, sys
from glob import glob
from Bparker.Bparker import return_B as Bparker_vector
from numpy.linalg import norm


def Lc_memilia(r=1.0):
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
	HYDRA cluster only.
    """
    fout = h5(fname_out, 'w')
    for r in range(wsize):
        fnm_inp = fname_out+'_%02d'%r
        finp = h5(fnm_inp, 'r')
        cont = finp.keys()       # list of groups
        for c in cont:           # iterate over each group
            finp.copy(c, fout)
        finp.close()
	print " --> removing partial files..."
        os.system('rm {fname}'.format(fname=fnm_inp))

    for pnm in psim.keys():
        fout['psim/%s'%pnm] = psim[pnm]

    print " ----> We generated: "+fout.filename
    fout.close()
    print " ----> FINISHED UNIFYING OUTPUT :D"



class Hmgr:
    """ class to unify all histograms of 
    step-sizes of every simulation with 
    several particles and one B-realization """
    def __init__(self, file, nbin=1000):
        self.f = file
        bmin, bmax = self.get_hstep_extremes()
        self.dbin = dbin = (bmax-bmin)/nbin
        self.hbin = np.arange(bmin+0.5*dbin, bmax, dbin)
        self.nbin = nbin
        # build histogram bounds
        self.hbd = hbd = np.zeros(self.nbin+1)
        hbd[:-1] = self.hbin - 0.5*self.dbin
        hbd[-1] = self.hbin[-1] + 0.5*self.dbin
        # initialize histogram in zero counts
        self.h = np.zeros(nbin)

    def get_hstep_extremes(self):
        """ Obtain the max && min of all the histogram domains """
        PNAMES  = self.f.keys()
        self.Np = len(PNAMES)
        self.bmin = 1.0e31
        self.bmax = 0.0
        for pnm in PNAMES:
            if not pnm.startswith('pla'):
                continue
            print '--->', pnm
            hx = self.f[pnm+'/HistStep/bins_StepPart'].value
            self.bmin = min(self.bmin, hx[0])
            self.bmax = max(self.bmax, hx[-1])
        return self.bmin, self.bmax

    def pile_to_hist(self, hx, hin):
        hbd = self.hbd
        #assert (hx[0]>=hbd[0]) & (hx[-1]<=hbd[-1]), \
        #    " RIDICULOUS ***ERROR*** IN OUR BOUNDARIES!!"

        for i in range(hx.size):
            #i = 0
            cc = (hx[i]>hbd[:-1]) & (hx[i]<hbd[1:])
            indx = find(cc)
            self.h[indx] += hin[i]



def load_traj(fname):
    f      = h5(fname, 'r')
    print " ---> reading: " + fname
    PNAMES = f.keys()    # particle names
    n      = len(PNAMES) # nmbr of plas in this file
    # take a sample to find out the times
    tadim  = f[PNAMES[0]+'/tadim'].value
    nt     = tadim.size 
    x = np.zeros((n,nt))
    y = np.zeros((n,nt))
    z = np.zeros((n,nt))
    for pname, i in zip(PNAMES, range(n)):
        if not pname.startswith('pla'):
            continue

        print " --> pname: ", pname
        x[i,:], y[i,:], z[i,:] = f[pname+'/xyz'].value.T # [1]

    wc   = f[PNAMES[0]+'/scl_wc'].value  # [s-1]
    rl   = f[PNAMES[0]+'/scl_rl'].value  # [cm]
    beta = f[PNAMES[0]+'/scl_beta'].value  # [1]
    o = {
    'x': x*rl, 'y':y*rl, 'z': z*rl,  # [cm] (n, nt)
    'tadim': tadim,
    'wc': wc, 'rl': rl, 'beta': beta,
    'nplas': n, 'ntime': nt,
    }
    return o


def get_sqrs(fname_inp):
    o     = load_traj(fname_inp)
    nt    = o['ntime']
    nplas = o['nplas']
    x, y, z = o['x'], o['y'], o['z']  # [cm] (nplas, nt)
    rl    = o['rl']                   # [cm]
    wc    = o['wc']                   # [s-1]
    tadim = o['tadim']                # [1]
    tdim  = tadim/wc                  # [s]
    AUincm = 1.5e13                   # [cm]
    # promediamos sobre particulas
    x2 = (x*x).mean(axis=0)           # [cm^2] 
    y2 = (y*y).mean(axis=0)           # [cm^2] 
    z2 = (z*z).mean(axis=0)           # [cm^2] 
    kxx = x2/tdim                     # [cm2/s]
    kyy = y2/tdim                     # [cm2/s]
    kzz = z2/tdim                     # [cm2/s]
    out = {
    'kxx': kxx, 'kyy': kyy, 'kzz': kzz,
    'tdim': tdim, 'tadim':tadim, 
    'wc': wc, 'rl': rl, 
    'nplas': nplas, 'ntime': nt,
    }
    return out


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


#w2file(rank, 'test.h5', a, 'a%02d'%rank)
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



def SaveToFile(m, dpath='', f=None, nbin=None):
    """
    input:
    - m     : cython-wrapper-module
    - dpath : directory-path for .h5 file
    - f     : output file object
    - nbin  : number of bins for step-size histos
    """
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
    #--- variables
    for name in vsave.keys():
        dname = dpath + name   # anteponemos algo
        f[dname] = vsave[name]

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
    f[dpath+'HistStep/HStep'] = HStp.sum(axis=0)
    f[dpath+'HistStep/bins_StepPart'] = bins_StepPart
    f[dpath+'HistStep/nbin'] = nbin



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
