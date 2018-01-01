#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
""" decodes an hex to a hash-key, so that I can
deduce the associated .h5 files, and list them.
"""

import shared.funcs_post as sfp
import os, sys
import argparse, h5py, logging

#--- retrieve args
parser = argparse.ArgumentParser(
description="""
Decodes an hex to a hash-key, so that I can
deduce the associated .h5 files, and list them.
IMPORTANT: the generated hash key (associated to the list of
input files) is only based on the cardinal numbers of the .h5 files; 
but NOT in the whole basenames.
""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-hash', '--hash',
type=str,
#default='<hash>'.format(**os.environ),
help='hash code that identifies the set of .h5 analyzed files. Enter only the hash (without the ".key" string)',
)
parser.add_argument(
'-content', '--list_contents',
#type='store_',
action='store_true',
default=False,
help='whether or not to list parameter values of the HDF5 files',
)
parser.add_argument(
'-dirkey', '--dir_key',
type=str,
default='<PLA1>/figs' if 'PLA1' not in os.environ else '{PLA1}/figs'.format(**os.environ),
help='directory where the hash-key files (*.key) are stored',
)
parser.add_argument(
'-dirh5', '--dir_h5',
type=str,
default='<PLA1>/out' if 'PLA1' not in os.environ else '{PLA1}/out'.format(**os.environ),
help='directory where the input .h5 files are stored',
)

pa = parser.parse_args()

#hexsuffix = sys.argv[1]
hexsuffix = pa.hash

#dir_key   = pa.dir_key #'{PLA1}/figs'.format(**os.environ)
fname_inp = '{dir}/{hsuff}.key'.format(dir=pa.dir_key, hsuff=hexsuffix)
if not os.path.isfile(fname_inp):
    " -----> no existe el archivo .key!" 
    exit()

#dir_src = pa.dir_h5 #'{PLA1}/out'.format(**os.environ)
decoded_str, IDs, prefix = sfp.DecodeHex_and_GetIDs(fname_inp)

ok = True
for ID in IDs:
    # we assume file IDs have 3 digits (see below: %03d)
    fnm = pa.dir_h5 + '/' + prefix + '%03d'%ID + '.h5'
    os.system('ls -lhtr '+fnm)
    ok &= os.path.isfile(fnm)

if ok:
    print " ---> all .h5 files exist!"
else:
    print " ---> *not* all files exist!"

if ok and pa.list_contents:
    from tabulate import tabulate
    p_comm, p_diff = sfp.compare_psim(pa.dir_h5, prefix, IDs, loglevel=logging.WARNING)

    #--- primero, hacemos la tabla para los parametros en comun
    tbcomm = [[ 'name', 'value' ]]
    # primero las semillas
    for nm in p_comm.keys():
        if not nm.startswith('sem_'):
            tbcomm += [[ nm, p_comm[nm] ]]

    tab_comm = tabulate(tbcomm, tablefmt='simple', headers='firstrow')
    print("\n ===== PARAMETROS COMUNES ===== ")
    print(tab_comm)

    #--- ahora hacemos la tex-tabla para los parametros diferentes
    header = []
    tbdiff = []
    for _myid in IDs:
        header += [ '%03d'%_myid ]

    tbdiff += [ header ]
    for nm in p_diff.keys():
        pars = []
        for myid in IDs:
            fname_inp = pa.dir_h5 + '/' + prefix + '%03d.h5' % myid
            with h5py.File(fname_inp, 'r') as f:
                par = f['psim/'+nm].value
                pars += [ '%2.2e' % par ]
        #tbdiff += [ ['\\texttt{$%s$}'%nm.replace('_','\_'),] + pars ]
        tbdiff += [ [nm,] + pars ]
    
    tab_diff = tabulate(tbdiff, tablefmt='simple', headers='firstrow')
    print("\n ===== PARAMETROS DIFERENTES ===== ")
    print(tab_diff)
    


#EOF
