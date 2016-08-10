#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
What we do here:
- process .h5 output-files to generate .key (hex of hash generated 
code) && .pdf (figures).
- example:
$ python ./juguete.py --IDs 37,38,39,40,41 --legend Nm_slab,Nm_2d
"""
import os, sys
from os.path import isfile, isdir
import shared.funcs_post as sfp
import argparse

# retrieve the IDentifiers :-)
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--IDs', type=str)
parser.add_argument('-l', '--legend', type=str)
parser.add_argument('-p', '--prefix', type=str, default='o_')
try:
    pa = parser.parse_args()
    ids = pa.IDs
    mylist = map(int, ids.split(','))
    mylabels = pa.legend.split(',')
    prefix = pa.prefix
    print pa
except IOError, msg:
    parser.error(str(msg))

ps = {
'dir_src'   : '%s/out' % os.environ['PLA1'],
'dir_dst'   : '%s/figs' % os.environ['PLA1'],
'id'        : mylist,   # list of files identifiers
'label'     : mylabels, # list of params for legend
}

ga = sfp.GenAnalysis(ps, prefix=prefix)
ga.gen_hash()  # genera el identificador
ga.make_pdf()

#EOF
