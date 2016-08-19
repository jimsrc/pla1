#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
""" decodes an hex to a hash-key, so that I can
deduce the associated .h5 files, and list them.
"""

import shared.funcs_post as sfp
import os, sys

hexsuffix = sys.argv[1]

dir_src   = '{PLA1}/figs'.format(**os.environ)
fname_inp = '{dir}/{hsuff}.key'.format(dir=dir_src, hsuff=hexsuffix)
if not os.path.isfile(fname_inp):
    " -----> no existe el archivo .key!" 
    exit()

dir_src = '{PLA1}/out'.format(**os.environ)
decoded_str, IDs, prefix = sfp.DecodeHex_and_GetIDs(fname_inp)

ok = True
for ID in IDs:
    # we assume file IDs have 3 digits (see below: %03d)
    fnm = dir_src + '/' + prefix + '%03d'%ID + '.h5'
    os.system('ls -lhtr '+fnm)
    ok &= os.path.isfile(fnm)

if ok:
    print " ---> all .h5 files exist!"
else:
    print " ---> *not* all files exist!"




##--- decodification method
#nID = len(decoded)/4
#ok = True
#for i in range(nID):
#    #TODO: don't always assume 'o_'
#    fnm = '{PLA1}/out'.format(**os.environ) + '/o_' + decoded[4*i:4*(i+1)]+'.h5'
#    os.system('ls -lhtr '+fnm)
#    ok &= os.path.isfile(fnm)
#
#if ok:
#    print " ---> all .h5 files exist!"
#else:
#    print " ---> *not* all files exist!"
#EOF
