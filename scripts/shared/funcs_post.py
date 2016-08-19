#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from h5py import File as h5
from os.path import isfile, isdir
import numpy as np
from pylab import (
    pause, find, figure, close
)
#from numpy import min, max
import os, sys
from glob import glob
from Bparker.Bparker import return_B as Bparker_vector
from numpy.linalg import norm
# for hash generation/encoding
from Crypto.Cipher import AES
import base64
# for multiple-page pdf generation
from matplotlib.backends.backend_pdf import PdfPages
# para hacer tablas, y transformarlas en .tex
from tabulate import tabulate
# para mergear .pdf files
from PyPDF2 import PdfFileMerger, PdfFileReader
from funcs import Bo_parker, Lc_memilia

#--- global constants
M_PI = np.pi
M_E  = np.e
AUincm = 1.5e13                   # [cm]


def DecodeHex_and_GetIDs(fname_key=None):
    """
    Get the list of id-file-numbers from a 
    string stored inside the file 'fname_key'.
    IMPORTANT: this should be upated with the
    codification method implemented in
    the GenAnalysis::gen_codification() routine.
    """
    f = open(fname_key, 'r')
    key = f.readline()
    prefix = f.readline()[:-1] # avoiding the '\n'
    gh = GenHash()
    decoded_str = gh.decode(encoded=key)
    #--- here is the decoding method
    #n = len(decoded_str)/4
    #IDs = []
    #for i in range(n):
    #    IDs += [ int(decoded_str[4*i:4*(i+1)]) ]
    M = int(decoded_str)
    m = 10 # NOTE: this MUST BE the same as 
           # in GenAnalysis::gen_codification() !!
    n, IDs, flag = 0, [], True
    while flag:
        an = int( M/(2**(n*m)) )%(2**m) 
        if an==0:
            flag = False
        else:
            IDs += [ an ]
            n += 1
    #-------------------------------
    return decoded_str, IDs, prefix


class GenHash(object):
    # the block size for the cipher object; must be 16, 24, or 32 for AES
    BLOCK_SIZE = 32
    # padding character, so that the value encrypted is 
    # always multiple of 'BLOCK_SIZE'
    PADDING = '{'
    # one-liner to sufficiently pad the text to be encrypted
    pad = lambda self, s:\
        s + (self.BLOCK_SIZE - len(s) % self.BLOCK_SIZE) * self.PADDING
    # one-liners to encrypt/encode and decrypt/decode a string
    # encrypt with AES, encode with base64
    EncodeAES = lambda self, c, s:\
        base64.b64encode(c.encrypt(self.pad(s)))
    DecodeAES = lambda self, c, e:\
        c.decrypt(base64.b64decode(e)).rstrip(self.PADDING)

    def __init__(self, **kargs):
        """
        input:
        - secret: used for encryption
        routine:
        - create a "secret" key (used to encrypt the message/key)
        - create a cipher object (using the "secret")
        """
        # generate a random secret key (if not given as argument)
        #secret = kargs.get('secret', os.urandom(self.BLOCK_SIZE))
        secret = '___JimmyMasias__' # THIS MUST BE FIXED!!
        self.myhash = {}
        # create a cipher object using the random secret
        self.myhash['cipher'] = AES.new(secret)

    def encode(self, mykey):
        """ encode a string """
        self.myhash['key'] =  mykey

        self.myhash['encoded'] = self.EncodeAES(self.myhash['cipher'], mykey)
        print 'Encrypted string:\n', self.myhash['encoded']
        return self.myhash['encoded']

    def decode(self, **kargs):
        """ decode the encoded string """
        encoded = kargs.get('encoded', self.myhash.get('encoded',None))
        assert encoded is not None,\
            " ### ERROR ###: Need 'encode' string from somewhere!"
        decoded = self.DecodeAES(self.myhash['cipher'], encoded)
        print 'Decrypted string:\n', decoded
        return decoded


class GenAnalysis(object):
    def __init__(self, ps, **kargs):
        self.idlist = ps['id']
        self.ps     = ps
        #---- other arguments
        self.fprefix = kargs.get('prefix', 'o_') # input fname prefix

    def gen_codification(self):
        """
        Generates a string code based on our
        list 'self.idlist' of id-numbers
        Returns a string.
        """
        # encode a string
        m  = 10 # (*)
        nf = len(self.idlist)
        M  = 0
        for n in range(nf):
            M += self.idlist[n]*2**(n*m)
        MyKey = str(M)
        # (*) so we can store list of numbers whose 
        #     values doesn't exceed 2**m, for each of them.
        return MyKey

    def gen_hash(self):
        gh = GenHash()
        # I want hash properties available
        self.myhash = gh.myhash
        # build string codification
        MyKey = self.gen_codification() #'This_passWD :P'
        return gh.encode(MyKey)

    def make_pdf(self):
        """ 
        We need:
        - hash identifier 'self.myhash["encoded"]' (e.g. see 
        the self.gen_hash() routine).
        What we do:
        Llama a los generadores de figuras 'GralPlot::plot_...()', y los
        pone uno en c/pagina de .pdf
        Antes de esto, genera una pagina q contiene la tabla de 
        simulation-parameters, y la pone antes de los plots.
        """
        ps = self.ps
        fname_base = self.myhash['encoded'].encode('hex')[:16]

        # construimos la pagina de parametros! :-)
        PdfOk, fname_param = self.build_params_pdf(fname_base)
        print " ---> PdfOk: ", PdfOk 
        
        fname_out_tmp  = '_tmp_fig_' + fname_base + '.pdf'
        fname_out = self.ps['dir_dst'] + '/fig_' + fname_base + '.pdf'
        pdf_pages = PdfPages(fname_out_tmp)

        gp = GralPlot(ps, prefix=self.fprefix, check=None, check_all=None)
        gp.do_checks()
        gp.build_labels()

        #--- 1st page
        fig, ax = gp.plot_errVel(OneFigFile=False)
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- 2nd page
        fig, ax, ax2 = gp.plot_errdy(OneFigFile=False)
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- 3rd page
        Ks = ('xx', 'yy', 'zz')
        for kk in Ks: #--- plot kxx, kyy, kzz
            fig, ax = gp.plot_kdiff(xaxis='mfp', nm=kk)
            pdf_pages.savefig(fig, bbox_inches='tight')
            close(fig)

        #--- 4th page
        fig, ax = gp.plot_TauColl()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- 5th page
        fig, ax = gp.plot_HistThetaColl()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- Write the PDF document to the disk
        pdf_pages.close()
        print " ---> we generated the temp-file: " + fname_out_tmp

        # ahora mergueamos los .pdf
        fnames_to_mergue = [fname_param, fname_out_tmp]
        merger = PdfFileMerger()
        for fname in fnames_to_mergue:
            merger.append(PdfFileReader(file(fname, 'rb')))

        merger.write(fname_out)
        print " -----> Mergueamos los .pdf aqui:\n " + fname_out 
        
        # borrando cosas temporales
        fnames_gb = fnames_to_mergue + [self.fname_tab_base+'*']
        print " -----> eliminandos .pdf temporales: ", fnames_gb
        for fnm in fnames_gb:
            os.system('rm '+fnm)

        #--- save code into an ASCII .key file (with my identifier)
        fname_key = self.ps['dir_dst'] + '/' + fname_base + '.key'
        os.system('echo ' + self.myhash['encoded'] + ' > '+fname_key)
        os.system('echo ' + self.fprefix + ' >> '+fname_key)
        print " ---> saved key (and prefix) into:\n"+fname_key

    def build_params_pdf(self, fname_base):
        p_comm, p_diff = self.compare_psim() # dictionaries
        tbcomm = [['name', 'value']]
        tbdiff = []

        #--- primero, hacemos la tabla para los parametros en comun
        # primero las semillas
        for nm in p_comm.keys():
            if nm.startswith('sem_'):
                name = '\\texttt{%s}'%nm.replace('_', u'\_')
                tbcomm += [[ name, p_comm[nm] ]]
        # ahora si el resto
        for nm in p_comm.keys():
            if not nm.startswith('sem_'):
                name = '\\texttt{%s}'%nm.replace('_', u'\_')
                tbcomm += [[ name, p_comm[nm] ]]

        # build tex table
        TexTab_comm = tabulate(tbcomm, tablefmt='latex', headers='firstrow')

        #--- ahora hacemos la tex-tabla para los parametros diferentes
        header = [u'name \\texttt{\\textbackslash} ID',]
        for myid in self.idlist:
            header += [ '%04d'%myid ]

        tbdiff += [ header ]
        for nm in p_diff.keys():
            pars = []
            for myid in self.idlist:
                fname_inp = self.ps['dir_src'] + '/' + self.fprefix + '%04d.h5' % myid
                with h5(fname_inp, 'r') as f:
                    par = f['psim/'+nm].value
                    pars += [ '%2.2e' % par ]
            tbdiff += [ ['\\texttt{$%s$}'%nm.replace('_','\_'),] + pars ]
        
        TexTab_diff = tabulate(tbdiff, tablefmt='latex', headers='firstrow')

        #--- beginin of .tex document
        with open(os.environ['HOME']+'/utils/tex_begin.txt', 'r') as f:
            tex_begin_lines = f.readlines()

        self.fname_tab_base = fname_tab_base = '_tmp.table_' + fname_base
        f = open(fname_tab_base+'.tex', 'w')
        for line in tex_begin_lines:
            f.write(line)
        # table of common params
        f.write(u'\\\\ \n {\\bf Common sim-parameters} \\\\ \n')
        for line in TexTab_comm:
            f.write(line)
        # table of different params
        f.write('\\vspace{1cm} \n')
        f.write(u'\\\\ \n {\\bf Different sim-parameters} \\\\ \n')
        for line in TexTab_diff:
            f.write(line)

        f.write('\n\end{document}')
        f.close()
        #--- end of .tex document

        cmd = 'pdflatex --interaction=nonstopmode {fname}'.format(fname=fname_tab_base+'.tex')
        return os.system(cmd), fname_tab_base+'.pdf'

    def compare_psim(self):
        """ Compare simulation-parameters to identify which
        are the same in each file, and which are different.

        output:
        - p_comm    : parameters in common
        - p_diff    : parameters different from one file to another
        """

        p_comm = {} # dict de params iguales
        p_diff = {} # dict de parametros en q difieren
        
        fname_inp_h5_ = self.ps['dir_src'] + '/' + self.fprefix + '%04d.h5' % self.idlist[0] # pick one
        f = h5(fname_inp_h5_, 'r')
        p_test = {} # dict de prueba, para comparar si todos son iguales
        for pnm in f['psim'].keys():
            p_test[pnm] = [ f['psim/'+pnm].value ] 

        for myid in self.idlist:
            fname_inp_h5 = self.ps['dir_src'] + '/' + self.fprefix + '%04d.h5' % myid
            f = h5(fname_inp_h5, 'r')
            print " ------- " + fname_inp_h5 + " ------- "
            for pnm in f['psim'].keys():
                if pnm in ('th', 'mu'):
                    continue

                fpar = f['psim/'+pnm].value # parameter from file
                if fpar == p_test[pnm]:
                    p_comm[pnm] = fpar
                else:
                    p_diff[pnm] = fpar

        print " ############ pars in common:"
        for nm in p_comm.keys():
            print nm, p_comm[nm]
        print " ############ pars different:"
        for nm in p_diff.keys():
            print nm, p_diff[nm]

        self.pars = {}
        self.pars['common'] = p_comm
        self.pars['different'] = p_diff

        return p_comm, p_diff


class GralPlot(object):
    def __init__(self, ps, prefix='o_', check=None, check_all=False):
        self.ps = ps
        self.prefix = prefix
        self.check = check
        self.check_all = check_all
        # symbols to iterate over
        self.sym = ('o', 's', '^', '*')

    def build_labels(self):
        """ construye los labels (legenda) de acuerdo a
        la lista de parametros en 'self.ps['label'].keys()'
        """
        #--- nombres/labels de los parametros
        names = ''
        if len(self.ps['label'])==1:
            names = self.ps['label'][0] + ': '
        else:
            for name in self.ps['label']:
                names += name + ', '
            names += ': '
        
        #--- valores de los parametros 
        values = {}
        for fid in self.ps['id']:
            fname_inp = self.ps['dir_src']+'/'+self.prefix+'%04d'%fid+'.h5'
            f = h5(fname_inp, 'r')
            values[fid] = ''
            for name in self.ps['label']:
                values[fid] += '$%G$' % f['psim/'+name].value
                if len(self.ps['label'])>1:  # si hay mas de uno,
                    values[fid] += ', '      # pongo coma.

        self.MyLabels = {}
        for fid in self.ps['id']:
            self.MyLabels[fid] = names + values[fid]

    def do_checks(self):
        """ reviso q ciertos keys existan en el output .h5 """
        ps = self.ps
        # should we check file keys?
        if self.check is not None:
            for fid in ps['id']:
                fname_inp = ps['dir_src']+'/'+self.prefix+'%04d'%fid+'.h5'
                f = h5(fname_inp, 'r')
                for ch in self.check:
                    # check we have that key in file
                    ok = ch in f['psim'].keys()
                    assert ok, ' ---> MISSING KEY in .h5 file!!'

        # should we check *all* 'ps' keys?
        if self.check_all:
            for fid, i in zip(ps['id'], range(len(ps['id']))):
                fname_inp = ps['dir_src']+'/'+self.prefix+'%04d'%fid+'.h5'
                f = h5(fname_inp, 'r')
                for ch in ps.keys(): # check all parameters
                    if ch.startswith(('dir_','label','id')): # skip
                        continue
                    # check we have that key in file
                    ok = ch in f['psim'].keys()
                    if not ok:
                        print "---> MISSING KEY in .h5 file: ", ch
                        raise SystemExit
                    #print " -----> i: ", i , ch
                    if ps[ch][i] != f['psim/'+ch].value:
                        print " ---> different simul-param: ", ps[ch][i], f['psim/'+ch].value, ch, fid
                        print " ---> Aborting..."
                        raise SystemExit
            f.close()
            print "\n ####################################### "
            print "      ALL SIM-PARAMETERS CHECK OK!"
            print " #######################################\n"

    def plot_errVel(self, OneFigFile=False, xlim=None, ylim=None):
        """
        - ylim: tuple for ax.set_ylim()
        """
        ps  = self.ps
        sym = self.sym
        ##--- figure
        fig = figure(1, figsize=(6,4))
        ax  = fig.add_subplot(111)
        #--- figname
        FigCode = ''
        for myid in ps['id']: FigCode += '%04d'%myid
        fname_fig = ps['dir_dst'] + '/errEk_' + FigCode + '.png'

        # iterate over all input-files
        id_indexes = range(len(ps['id'])) # indexes
        err_min, err_max = 1e31, -1e31 # para ajustar el set_ylim()
        for fid, i in zip(ps['id'], id_indexes):
            fname_inp = ps['dir_src']+'/'+self.prefix+'%04d'%fid+'.h5'
            f = h5(fname_inp, 'r')
            tadim = f['pla000/tadim']
            PNAMES = [] # particle names
            for nm in f.keys():
                if nm.startswith('pla'):
                    PNAMES += [ nm ]

            Np = len(PNAMES)
            err = np.zeros((Np,tadim.size))
            for pnm, ip in zip(PNAMES, range(Np)):
                err[ip,:] = np.abs(f[pnm+'/err'].value)

            #--- stats over pla realizations
            err_med = np.median(err, axis=0)
            err_avr = err.mean(axis=0)
            err_std = err.std(axis=0)
            #--- calcula los minimos sin tomar en cuenta el 1er tiempo
            err_min = np.min([err_min, np.min(err_avr[1:])])
            err_max = np.max([err_max, np.max(err_avr[1:])])

            isym = np.mod(i,len(sym))
            opt = {'ms': 3, 'mec':'none', 'marker': sym[isym-1],'ls':''}
            label = self.MyLabels[fid]  #ps['label'][i]
            ax.plot(tadim, err_avr, label=label, **opt)

        ax.legend(loc='best', fontsize=7)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.grid(True)
        ax.set_ylabel('energy error  [1]')
        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim(err_min, err_max)

        if OneFigFile:
            fig.savefig(fname_fig, dpi=200, bbox_inches='tight')
            close(fig)
        else:
            return fig, ax

    def plot_errdy(self, OneFigFile=False, nbin=1000, **kargs):
        """
        input:
        - nbin: number of bins for step-size histogram
        """
        AUincm = 1.5e13             # [cm]
        ps  = self.ps
        sym = self.sym
        xlim = kargs.get('xlim', None)
        ylim = kargs.get('ylim', None)
        #--- figure
        fig = figure(1, figsize=(6,4))
        ax  = fig.add_subplot(111)
        ax2 = ax.twiny()

        #--- figname
        FigCode = ''
        for myid in ps['id']: FigCode += '%04d'%myid
        fname_fig = ps['dir_dst'] + '/errdy_' + FigCode + '.png'

        # iterate over all input-files
        id_indexes = range(len(ps['id'])) # indexes
        for fid, i in zip(ps['id'], id_indexes):
            fname_inp = ps['dir_src']+'/'+self.prefix+'%04d'%fid+'.h5'
            f = h5(fname_inp, 'r')
            PNAMES = f.keys()

            hmg = Hmgr(f, nbin=nbin)
            hmg.get_hstep_extremes()

            for pnm in PNAMES:
                if not pnm.startswith('pla'):
                    continue

                h  = f[pnm+'/HistStep/HStep'].value
                hx = f[pnm+'/HistStep/bins_StepPart'].value
                hmg.pile_to_hist(hx, h)

            isym = np.mod(i,len(sym))
            opt = {'ms': 3, 'mec':'none', 'marker': sym[isym-1], 'ls':''}
            label = self.MyLabels[fid] #ps['label'][i]
            lmin = f['psim/lmin'].value # [r_larmor]
            dRbin = hmg.hbin/lmin # (eq. np-1)
            ax2.plot(dRbin, hmg.h, label=label, **opt)
            ax2.set_xlim(dRbin[0], dRbin[-1])
            ax2.set_xlabel('$\Delta r/\lambda_{min}$')
            ax2.set_xscale('log')
            ax.plot(hmg.hbin, hmg.h, label=label, **opt)
            ax.legend(loc='best', fontsize=7)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel('$\Omega dt$')
            ax.set_ylabel('#')
            ax.set_xlim(hmg.hbin[0], hmg.hbin[-1])
            if ylim is not None:
                ax.set_ylim(ylim)

            f.close()

        ax.grid(True)
        if OneFigFile:
            fig.savefig(fname_fig, dpi=200, bbox_inches='tight')
            close(fig)
        else:
            return fig, ax, ax2

    def plot_kdiff(self, xaxis='mfp', nm='xx', **kargs):
        assert xaxis in ('kdiff','mfp'),\
            ' ---> ERROR: wrong xaxis:%s @plot_kdiff()'%xaxis
        kk     = 'k'+nm # kxx,kyy,or kzz
        ps     = self.ps
        o      = {}
        sym    = self.sym
        #--- extra args
        xlim   = kargs.get('xlim', None)
        ylim   = kargs.get('ylim', None)
        xscale = kargs.get('xscale', 'log')
        yscale = kargs.get('yscale', 'log')
        #--- plot kxx, kyy, kzz
        print " ---> plotting " + kk + ':'
        #--- figura
        fig = figure(1, figsize=(6,4))
        ax  = fig.add_subplot(111)
        # iterate over all input-files
        id_indexes = range(len(ps['id'])) # indexes
        for fid, i in zip(ps['id'], id_indexes):
            fname_inp = ps['dir_src']+'/'+self.prefix+'%04d'%fid+'.h5'
            o[kk] = get_sqrs(fname_inp) # w/ corrected dimensions
            with h5(fname_inp, 'r') as f:
                Lc_s = f['psim/Lc_slab'].value # [1]
            tadim = o[kk]['tadim']      # [1]  (correct dimensions)
            kprof = o[kk][kk]           # [1]  (")
            isym  = np.mod(i,len(sym))
            print " i, len(sym), isym: ", i, len(sym), isym
            opt = {
            'ms'        : 2,
            'lw'        : 0.5,
            'marker'    : sym[isym-1],
            'label'     : self.MyLabels[fid],
            'alpha'     : 0.6,
            'mec'       : 'none',
            }
            mfp = 3.*kprof/Lc_s # [1] normalized (correct) units
            ax.plot(tadim, mfp, '-o', **opt)

        ax.set_yscale(yscale)
        ax.set_xscale(xscale)
        ax.set_xlabel('$\Omega t$ [1]')
        ax.set_ylabel('$\lambda_{%s}/L_c^{slab}$ [1]'%nm)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.legend(loc='best', fontsize=6)
        ax.grid()
        return fig, ax # dicts for (kxx,kyy,kzz)

    def plot_HistThetaColl(self, **kargs):
        ps     = self.ps
        o      = {}
        sym = self.sym
        #--- extra args
        xlim   = kargs.get('xlim', None)
        ylim   = kargs.get('ylim', None)
        xscale = kargs.get('xscale', 'log')
        yscale = kargs.get('yscale', 'log')
        #--- figura
        fig = figure(1, figsize=(6,4))
        ax  = fig.add_subplot(111)
        # iterate over all input-files
        id_indexes = range(len(ps['id'])) # indexes
        for fid, i in zip(ps['id'], id_indexes):
            fname_inp = ps['dir_src'] + '/'+self.prefix+'%04d'%fid + '.h5'
            ht = HThetaColl(fname_inp)
            h = ht.SumHsts_over_plas() # my histograms!!
            if h is 0:
                ax.text(.5, .5, 'sorry, no histograms for step-sizes.',\
                  transform=ax.transAxes)
                continue # next 'fid'
            else:
                hx, hc = h['hbins'], h['hcnts']
            isym = np.mod(i,len(self.sym))
            opt = {'ms': 3, 'mec':'none', 'marker': sym[isym-1],'ls':''}
            label = self.MyLabels[fid]  #ps['label'][i]
            ax.plot(hx, hc, label=label, **opt)
        ax.legend(loc='best', fontsize=7)
        ax.set_yscale('linear')
        ax.grid(True)
        ax.set_ylabel('#')
        ax.set_xlabel('$\\theta_{coll}$ [deg]')
        return fig, ax

    def plot_TauColl(self, OneFigFile=False, xlim=None, ylim=None, **kargs):
        """ plot histogram of collision-times.
        """
        ps  = self.ps
        sym = self.sym
        yscale = kargs.get('yscale', 'log')
        ##--- figure
        fig = figure(1, figsize=(6,4))
        ax  = fig.add_subplot(111)
        # iterate over all input-files
        id_indexes = range(len(ps['id'])) # indexes
        err_min, err_max = 1e31, -1e31 # para ajustar el set_ylim()
        for fid, i in zip(ps['id'], id_indexes):
            fname_inp = ps['dir_src'] + '/'+self.prefix+'%04d'%fid + '.h5'
            f = h5(fname_inp, 'r')
            ht = HTauColl(fname_inp, nbin=1000)
            if ht is 0:
                ax.text(.5, .5, 'sorry, no histograms for collision-tau.',\
                transform=ax.transAxes)
                continue # next 'fid'

            h = ht.SumHsts_over_plas()
            hx, hc = h['hbins'], h['hcnts']
            # now lets plot :)
            isym = np.mod(i,len(self.sym))
            opt = {'ms': 3, 'mec':'none', 'marker': sym[isym-1],'ls':''}
            label = self.MyLabels[fid]
            hx_ = hx - np.log10(2.*M_PI) #hx_: log10(omega*tau/2pi)
            ax.plot(hx_, hc, label=label, **opt)
        ax.legend(loc='best', fontsize=7)
        ax.set_yscale(yscale)
        ax.grid(True)
        ax.set_ylabel('#')
        ax.set_xlabel('$log_{10}(\Omega \\tau_{coll})$')
        return fig, ax


class HTauColl(object):
    def __init__(self, fname_inp, nbin=1000, **kargs):
        self.fname_inp = fname_inp
        self.nbin = nbin
        # main check
        nameHsts = 'HistTau_log' # key-name of data
        with h5(fname_inp, 'r') as f:
            if nameHsts not in f['pla000'].keys():
                return 0 # finish!
        # hallar los minimos de los dominios
        # de todos los histogramas
        hmin, hmax = self.get_hist_extremes()
        self.dbin = dbin = (hmax-hmin)/nbin
        self.hbin = np.arange(hmin+0.5*dbin, hmax, dbin)
        self.nbin = nbin
        # build histogram bounds
        self.bd = bd = np.zeros(self.nbin+1)
        bd[:-1] = self.hbin - 0.5*self.dbin
        bd[-1] = self.hbin[-1] + 0.5*self.dbin

    def get_hist_extremes(self):
        """ Obtain the max && min of all the 
        histogram domains """
        print " ---> getting hist-log(tau)-extremes..."
        f = h5(self.fname_inp, 'r')
        PNAMES  = f.keys()
        self.Np = len(PNAMES)
        self.bmin = 1.0e31
        self.bmax = 0.0
        nameHsts = 'HistTau_log' # key-name of data
        for pnm in PNAMES:
            if not pnm.startswith('pla'):
                continue
            hc_, hx = f[pnm+'/'+nameHsts].value
            self.bmin = min(self.bmin, hx[0])
            self.bmax = max(self.bmax, hx[-1])
        return self.bmin, self.bmax

    def SumHsts_over_plas(self):
        """ unifica/suma los histogramas de todas
        las particulas, para una B-realization"""
        f = h5(self.fname_inp, 'r')
        nameHsts = 'HistTau_log' # key-name of data
        # initialize histogram in zero counts
        hcnts = np.zeros(self.nbin, dtype=np.float32)
        bd    = self.bd
        for dpath in f.keys():
            if not dpath.startswith('pla'):
                continue
            hc, hx = f[dpath+'/'+nameHsts].value
            for i in range(hx.size):
                cc = (hx[i]>bd[:-1]) & (hx[i]<bd[1:])
                ix = find(cc)
                hcnts[ix] += hc[i]
        self.myhist = {
        'hbins' : self.hbin,
        'hcnts' : hcnts,
        }
        return self.myhist



class HThetaColl(object):
    def __init__(self, fname_inp, nbin=1000, **kargs):
        self.fname_inp = fname_inp
        """ hacer algo q me diga en que
        intervalos de tiempo deberia armar
        los histogramas
        self.every = ..."""
        # main check
        nameHsts = 'HistThetaColl' # key-name of data
        with h5(fname_inp, 'r') as f:
            if nameHsts not in f['pla000'].keys():
                return 0 # finish!
        # hallar los minimos de los dominios
        # de todos los histogramas
        hmin, hmax = self.get_hist_extremes()
        self.dbin = dbin = (hmax-hmin)/nbin
        self.hbin = np.arange(hmin+0.5*dbin, hmax, dbin)
        self.nbin = nbin
        # build histogram bounds
        self.bd = bd = np.zeros(self.nbin+1)
        bd[:-1] = self.hbin - 0.5*self.dbin
        bd[-1] = self.hbin[-1] + 0.5*self.dbin

    def get_hist_extremes(self):
        """ Obtain the max && min of all the 
        histogram domains """
        print " ---> getting hist-theta-extremes..."
        f = h5(self.fname_inp, 'r')
        PNAMES  = f.keys()
        self.Np = len(PNAMES)
        self.bmin = 1.0e31
        self.bmax = 0.0
        nameHsts = 'HistThetaColl' # key-name of data
        for pnm in PNAMES:
            if not pnm.startswith('pla'):
                continue
            hc_, hx = f[pnm+'/'+nameHsts].value
            self.bmin = min(self.bmin, hx[0])
            self.bmax = max(self.bmax, hx[-1])
        return self.bmin, self.bmax

    def SumHsts_over_plas(self):
        """ unifica/suma los histogramas de todas
        las particulas, para una B-realization"""
        f = h5(self.fname_inp, 'r')
        nameHsts = 'HistThetaColl' # key-name of data
        # initialize histogram in zero counts
        hcnts = np.zeros(self.nbin, dtype=np.float32)
        bd    = self.bd
        for dpath in f.keys():
            if not dpath.startswith('pla'):
                continue

            hc, hx = f[dpath+'/'+nameHsts].value
            for i in range(hx.size):
                cc = (hx[i]>bd[:-1]) & (hx[i]<bd[1:])
                ix = find(cc)
                hcnts[ix] += hc[i]

        self.hist_stepsize = {
        'hbins' : self.hbin,
        'hcnts' : hcnts,
        }
        return self.hist_stepsize



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
        print " ---> getting hstep-extremes... "
        PNAMES  = self.f.keys()
        self.Np = len(PNAMES)
        self.bmin = 1.0e31
        self.bmax = 0.0
        for pnm in PNAMES:
            if not pnm.startswith('pla'):
                continue
            #print '--->', pnm
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
    # count plas
    PNAMES = []
    for pnm in f.keys():
        if pnm.startswith('pla'): PNAMES += [ pnm ]
    #PNAMES = f.keys()    # particle names
    n      = len(PNAMES) # nmbr of plas in this file
    # take a sample to know the times
    tadim  = f[PNAMES[0]+'/tadim'].value
    nt     = tadim.size 
    x = np.zeros((n,nt))
    y = np.zeros((n,nt))
    z = np.zeros((n,nt))
    for pname, i in zip(PNAMES, range(n)):
        x[i,:], y[i,:], z[i,:] = f[pname+'/xyz'].value.T # [1]

    o = {
    'x': x, 'y':y, 'z': z,  # [1] (n, nt)
    'tadim': tadim,
    'nplas': n, 'ntime': nt,
    }
    return o


def get_sqrs(fname_inp):
    o     = load_traj(fname_inp)
    nt    = o['ntime']
    nplas = o['nplas']
    x, y, z = o['x'], o['y'], o['z']  # [1] (nplas, nt)
    tadim = o['tadim']                # [1]
    #tdim  = tadim/wc                  # [s]
    #AUincm = 1.5e13                   # [cm]
    # promediamos sobre particulas
    x2 = (x*x).mean(axis=0)           # [1] 
    y2 = (y*y).mean(axis=0)           # [1] 
    z2 = (z*z).mean(axis=0)           # [1] 
    kxx = x2/(2.*tadim)               # [1]
    kyy = y2/(2.*tadim)               # [1]
    kzz = z2/(2.*tadim)               # [1]
    out = {
    'kxx': kxx, 'kyy': kyy, 'kzz': kzz, # [1]
    'tadim':tadim, 
    #'wc': wc, 'rl': rl, 
    'nplas': nplas, 'ntime': nt,
    }
    return out

#EOF
