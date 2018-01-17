#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import figure, close
import h5py, argparse, sys, os


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
parser.add_argument(
'-dirfig', '--dirfig',  # the .h5 file is assumed to have a '%03d' format for the ID
type=str,
help='output for all figures',
)
parser.add_argument(
'-tauhist', '--tauhist',
action='store_true',
default=False,
)
parser.add_argument(
'-indiv', '--indiv',
type=float,
nargs=2,
default=None,
metavar=('ThresLower', 'ThresUpper'),
help='range',
)
parser.add_argument(
'-meanprof', '--meanprof',
type=float,
nargs=2,
default=None,
metavar=('ThresLower', 'ThresUpper'),
help='range',
)
pa = parser.parse_args()


#+++++++++++++++++++++++++++++++++++++++++++++++++
pi2 = 2.0*np.pi
f = h5py.File(pa.fname_inp, 'r')
# let's see what kind of scattering times were captured
Npla=len([nm for nm in f.keys() if nm.startswith('pla')])


#+++++++++++++++++++++++++++++++++++++++++++++++++
if pa.tauhist:
    #taus =  []
    htau_nbin   = 100
    htau_range  = [0. ,4.]
    htau_res    = (htau_range[1] - htau_range[0])/htau_nbin
    htau        = np.zeros(100, dtype=np.int64)
    tau_bins    = np.arange(htau_range[0]+0.5*htau_res, htau_range[1], step=htau_res)
    for ipla in range(Npla):
        _taus = f['pla%03d/trails/tau_b'%ipla][...]/pi2
        #taus += [ f['pla%03d/trails/tau_b'%ipla][...] ]
        htau += np.histogram(_taus, bins=htau_nbin, range=htau_range, normed=False)[0]


    fig = figure(1, figsize=(6,4))
    ax  = fig.add_subplot(111)
    ax.plot(tau_bins, htau, '-o')
    #fig.show()
    ax.grid(1)
    ax.set_yscale('log')
    fig.savefig(pa.dirfig+'/taus_hits.png', bbox_inches='tight')
    close(fig)



#+++++++++++++++++++++++++++++++++++++++++++++++++
if pa.indiv is not None:
    taumin, taumax = pa.indiv
    assert taumax > taumin, ' [-] error here..\n'
    subdirfig = pa.dirfig+'/muprofs__%.2f_%.2f' % (taumin,taumax)
    os.system('mkdir -p '+subdirfig)
    zoom = 80
    for ipla in range(Npla):
        _taus = f['pla%03d/trails/tau_b'%ipla][...]/pi2
        cc = (_taus>taumin) & (_taus<taumax)
        # number of profiles with a coll-time in the range of interest
        n_muprofs = cc.nonzero()[0].size
        # read trajectories of all scatterings of this pla
        tr = f['pla%03d/trails/trajecs'%ipla][cc,:,:]
        fig = figure(1, figsize=(1*zoom,1*zoom))

        # iterate over all the bounces
        for ireb in range(n_muprofs):
            t, mu, xyz  = tr[ireb,:,0], tr[ireb,:,1], tr[ireb,:,2:]
            #x, y, z = xyz.T
            print " [*] plotting scatter/pla: %d/%d" % (ireb,ipla)
            # layout for this profile
            ax  = fig.add_subplot(n_muprofs/4+1, 4, ireb+1)
            ax.plot((t[0]-t)/pi2, np.abs(mu), '-')

            ax.axhline(y=0.0, lw=3., c='k', alpha=0.6)

            ax.grid(1)
            ax.set_ylim(0., 1.)

        # save plot
        figname = subdirfig + '/pla%03d.png'%ipla
        fig.savefig(figname, bbox_inches='tight')
        close(fig)
   


#+++++++++++++++++++++++++++++++++++++++++++++++++
if pa.meanprof:
    taumin, taumax = pa.meanprof
    assert taumax > taumin, ' [-] error here..\n'
    subdirfig = pa.dirfig+'/mean__%.2f_%.2f' % (taumin,taumax)
    os.system('mkdir -p '+subdirfig)
    zoom = 80
    for ipla in range(Npla):
        _taus = f['pla%03d/trails/tau_b'%ipla][...]/pi2
        cc = (_taus>taumin) & (_taus<taumax)
        # number of profiles with a coll-time in the range of interest
        #n_muprofs = cc.nonzero()[0].size
        # read trajectories of all scatterings of this pla
        tr = f['pla%03d/trails/trajecs'%ipla][cc,:,:]

        # average over all the bounces
        abs_mu  = np.abs(tr[:,:,1])     # absolute value
        mu_mean = abs_mu.mean(axis=0)
        mu_std  = abs_mu.std(axis=0)

        # grab a sample time array
        t    = tr[0,:,0]
        time = (t[0] - t)/pi2

        # figure
        print " [*] plot ipla: %d" % ipla
        fig = figure(1, figsize=(6,4))
        ax = fig.add_subplot(111)
        ax.plot(time, mu_mean, '-')
        ax.grid(1)
        ax.set_ylim(0., 1.)

        fig.savefig(subdirfig+'/pla%03d.png'%ipla, bbox_inches='tight')
        close(fig)


#EOF
