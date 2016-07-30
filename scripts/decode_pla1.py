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

with open(fname_inp, 'r') as f:
    myhash = f.readline()

gh = sfp.GenHash()

#myhash = hexstr.decode('hex')
decoded = gh.decode(encoded=myhash)
nID = len(decoded)/4
ok = True
for i in range(nID):
    fnm = '{PLA1}/out'.format(**os.environ) + '/o_' + decoded[4*i:4*(i+1)]+'.h5'
    os.system('ls -lhtr '+fnm)
    ok &= os.path.isfile(fnm)

if ok:
    print " ---> all .h5 files exist!"
#EOF
