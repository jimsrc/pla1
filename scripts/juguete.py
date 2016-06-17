#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import os, sys
from os.path import isfile, isdir
import shared.funcs as sf

ps = {
'dir_src'   : '%s/out' % os.environ['PLA1'],
'dir_dst'   : '%s/figs' % os.environ['PLA1'],
'id'        : (47, 48, 51),
'Nm_2d'     : (256, 128, 64),
'Nm_slab'   : (128, 128, 128),
'eps_o'     : (1e-5, 1e-5, 1e-6),
'label'     : ('uno', 'dos', 'tres'),
}

#ff.plotear_kdiff(ps=ps, check=('eps_o',), check_all=False)
gp = sf.GralPlot(ps=ps, check=None, check_all=False)
gp.do_checks()
gp.plot_kdiff()

#EOF
