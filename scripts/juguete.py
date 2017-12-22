#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
What we do here:
- process .h5 output-files to generate .key (hex of hash generated 
code) && .pdf (figures).
- example:
$ ./juguete.py -- --IDs 5,6,7,8 --legend Nm_slab,Nm_2d --prefix h_
"""
import os, sys
from os.path import isfile, isdir
import shared.funcs_post as sfp
import argparse

# retrieve the IDentifiers :-)
parser = argparse.ArgumentParser(
description="""
Process .h5 output-files to generate a .key (hex of hash generated 
code) && .pdf (figures).
The convention, assumed here, is than the input .h5 files are symlinks (not necessarily) to the real files. The aim with the symlinks is to categorize the original .h5 files; this way, a single ('original') .h5 file can belong to several categories of an analysis.
""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-i', '--IDs',  # the .h5 file is assumed to have a '%03d' format for the ID
type=str,
help='list of identifier numbers.',
)
parser.add_argument(
'-l', '--legend', 
type=str,
help='list of parameters to show in legend of figures.',
)
parser.add_argument(
'-p', '--prefix', 
type=str, 
default='h_',
help='prefix string of input filenames.',
)
parser.add_argument(
'-src', '--dir_src', 
type=str, 
default='%s/out' % os.environ['PLA1'],
help='input directory',
)
parser.add_argument(
'-dst', '--dir_dst',
type=str,
default='%s/figs' % os.environ['PLA1'],
help='output directory',
)
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
'dir_src'   : pa.dir_src,
'dir_dst'   : pa.dir_dst,
'id'        : mylist,   # list of files identifiers
'label'     : mylabels, # list of params for legend
}

ga = sfp.GenAnalysis(ps, prefix=prefix)
ga.gen_hash()  # genera el identificador
ga.make_pdf()

#EOF
