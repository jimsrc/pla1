"""
Module parsers, for the different magnetic field models.
"""
import numpy as np
import src.cython_wrapper as cw
import shared.funcs as sf
import mparsers as mp
#from pylab import pause
#--- parameters
#from params import (
#    nB, pd, psim, pother, AU_in_cm
#)
from h5py import File as h5
from os.path import isfile, isdir
import os, sys, argparse, time
from glob import glob
from numpy import power, log10

#--- globals
AUincm = 1.5e13        # [cm]
r2d    = 180./np.pi    # rad to degrees

# gral info
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
#++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++
# dB models
#++++++++++++++++++++++++++++++++++++++++++++++++++++
class dBmodel__slab2d(object):
    """
    Slab + 2D turbulence model
    """
    def __init__(self):
        """
        First thing to do is build the parser
        """
        self.help = """
        Module to read/process the output of the 'build_pressure.corr' stage.
        """
        self.parser = parser = argparse.ArgumentParser(
        description="""this gral description...""", 
        add_help=False
        )
        parser.add_argument(
        '-ro',
        type=float,
        default=0,
        help='heliodistance in AU.',
        )
        parser.add_argument(
        '-sigB', '--sigmaBoB',
        type=float,
        default=0.3,
        help='energy of the magnetic fluctuations relative to the guide field, (dB/Bo)^2.',
        )
        parser.add_argument(
        '-fslab',
        type=float,
        default=0.2,
        help='fraction for the Slab energy (relative to the total fluctuation dB^2)',
        )
        parser.add_argument(
        '-NmS',
        type=int,
        default=128,
        help='number of slab modes',
        )
        parser.add_argument(
        '-Nm2d',
        type=int,
        default=128,
        help='number of 2D modes',
        )
        parser.add_argument(
        '-lmin',
        type=float,
        nargs=2,
        default=[0,0],
        metavar=('LMIN-SLAB', 'LMIN-2D'),
        help="""minimum wavelengths (LMIN-SLAB, LMIN-2D; in units of the Larmor
        radii) of the Fourier spectrum. If this is used, then the -ro option must 
        not be used.""",
        )
        parser.add_argument(
        '-lmax',
        type=float,
        nargs=2,
        default=[0,0],
        metavar=('LMAX-SLAB', 'LMAX-2D'),
        help="""maximum wavelengths (LMAX-SLAB, LMAX-2D; in units of the Larmor 
        radii) of the Fourier spectrum. If this is used, then the -ro option must 
        not be used."""
        )
        parser.add_argument(
        '-Lc',
        type=float,
        default=0,
        help="""Correlation scale (not correlation lengths) for Slab, in units 
        of Larmor radii, for the turbulence model. If this is used, then 
        the -ro option must not be used."""
        )
        parser.add_argument(
        '-xi',
        type=float,
        default=1.0,
        help='ratio Lc_slab/Lc_2d, where Lc_slab is the correlation scale.',
        )

    def run(self, pa, **kws):
        #--- consistency checks
        # either specify the helioradius to determine some turb parameters, or
        # specify the parameters themself.
        assert (pa.ro==0) or (pa.lmin==[0,0] and pa.lmax==[0,0] and pa.Lc==0), \
        '\n [-] ilegal use of arguments: use -ro or -lmin .. -lmax .. -Lc .., not both!\n'
        # specify a number of picth angles, or a value of a pitch angle
        assert (pa.Nth==0) or (pa.tho is None), \
        '\n [-] either the -Nth or -tho can only be used, not both!\n' 


        #--- obtain the correlation length for this heliodistance 'ro'
        # correlation LENGTH in [AU]
        lc = sf.Lc_memilia(r=pa.ro) \
            if pa.ro else None
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
        # NOTE; the expression below comes from the above comments
        Lc_slab = power(10., 1.026*log10(lc) + 0.2755) \
            if pa.ro else None
            
        
        Rl = cw.calc_Rlarmor(       # [AU] Larmor radii
            rigidity=1.69604E+09, #1.69604E+09,    # [V]
            Bo=sf.Bo_parker(r=pa.ro)    # [Gauss]
            )/AUincm if pa.ro else None
        #--- set B-turbulence model
        pd = kws['pd']
        pd.update({
        'Nm_slab'        : pa.NmS,
        'Nm_2d'          : pa.Nm2d,
        'lmin_s'         : 5e-5/Rl  if pa.ro else pa.lmin[0], #[Rl] 
        'lmax_s'         : pa.ro/Rl if pa.ro else pa.lmax[0],  #[Rl] 
        'lmin_2d'        : 5e-5/Rl  if pa.ro else pa.lmin[1], #[Rl] 
        'lmax_2d'        : pa.ro/Rl if pa.ro else pa.lmax[1],  #[Rl] 
        'Lc_slab'        : Lc_slab/Rl if pa.ro else pa.Lc,  # [Rl] in units of Larmor-radii
        'xi'             : pa.xi, # [1] xi=Lc_2d/Lc_slab 
        'sigma_Bo_ratio' : pa.sigmaBoB,  # [1] fluctuation energy
        'ratio_slab'     : pa.fslab,     # [1] (energy_slab)/(energy_total)
        })
        #--- corregimos input
        psim = kws['psim']
        psim['tmax']     = pa.tmax     # [1/omega]
        lmin             = np.min([pd['lmin_s'], pd['lmin_2d']]) # smallest turb scale
        psim['atol']     = lmin*pa.eps # [1]
        psim['rtol']     = 0.0 #1e-6
        pother                  = kws['pother']
        pother['str_timescale'] = pa.tscale
        pother['nsave']         = pa.nsave

        #--- orientations of the initial velocities
        mu  = kws['mu']
        ph  = kws['ph']
        
        #--- output
        po = {}
        po.update(psim)
        po.update(pd)
        # add some stuff
        po.update({
        'r'     : pa.ro,        # [AU] heliodistance
        'eps_o' : pa.eps,       # [1]  precision
        'lmin'  : lmin,         # [1]  minimum turb scale (*)
        'RloLc' : Rl/Lc_slab if pa.ro else 1./pa.Lc,   # [1]  (r_larmor)/(Lc_slab)
        })
        
        #--- append trails info if given
        if pa.taubd is not None:
            po.update({
            'taubd' : np.array(pa.taubd), # list of borders (of bands)
            })
        
        # (*) according to the definitions above of 'lmin_*' they are
        # adimensional but they're in units of Rl.
        
        odir = '/'.join(pa.fname_out.split('/')[:-1])
        assert isdir(odir), '\n [-] doesnt exist output directory:\n %s\n'%odir
        
        #--- call simulator
        m = cw.mgr()
        m.set_Bmodel(pdict=pd, nB=kws['nB'])

        # MPI stuff
        #rank, wsize = kws['rank'], kws['wsize']
        comm        = kws['comm'] #MPI.COMM_WORLD
        rank        = comm.Get_rank()   # proc rank
        wsize       = comm.Get_size()   # number of proc

        #--- memory check
        # NOTE: we estimate the "requested data" in terms
        # of those pieces of C++ code where there is a
        # "WATCH_MEMORY" comment (see .cc sources!).
        STOT, avail_ram = sf.memory_stat(
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
        
        pla_bd = sf.equi_bounds(0, mu.size-1, wsize) # bounds
        plas   = np.arange(pla_bd[rank], pla_bd[rank+1]) # plas for each proc
        
        if rank==0 and isfile(pa.fname_out):  # backup if already exists
            os.system('mv {fname} {fname}_'.format(fname=pa.fname_out))

        fname_out_tmp = pa.fname_out+'_%02d'%rank # output of each processor
        fo = h5(fname_out_tmp, 'w') 
        nbin = 1000 # for step-size histograms
        for npla in plas: #[25:]:
            #--- set particle id && direction
            pother['i']         = npla
            if pa.taubd is not None:
                pother['__tau_bd'] = pa.taubd   # list of band borders
            psim['mu']          = mu[npla]
            psim['ph']          = ph[npla]
            print " mu, ph: %g/%g" % (mu[npla], ph[npla])
        
            m.build(**pother)
        
            m.SetSim(**psim)
            m.RunSim()
        
            print " [r:%03d] (pla:%d of %d) finished!" % (rank, npla, plas[-1])
            dpath = 'pla%03d/' % npla
            print " [r:%03d] writing: %s" % (rank, fo.filename)
            sf.SaveToFile(m, dpath, fo, nbin_step=nbin, btrails=pa.trails)
            m.clean()

        print " [r:%03d] closing: %s" % (rank, fo.filename)
        fo.close()
        time.sleep(1)
        os.system('touch %s_finished'%fname_out_tmp)
        print " [r:%03d] I'm finished!" % rank
        
        if rank==0:
            #--- let's check if everyone has already finished
            n = 0
            while n<wsize:
                time.sleep(2) # wait before next check
                fs = glob(pa.fname_out+'_*_finished') #list of files
                n  = len(fs)
            
            #--- now proceed to unify
            sf.unify_all(pa.fname_out, psim=po, wsize=wsize)
            #--- now we can delete the "*_finished" files
            for fnm in fs:
                os.system('rm '+fnm)


#EOF