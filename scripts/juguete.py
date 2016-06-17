#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import os, sys
from os.path import isfile, isdir
import shared.funcs as sf

ps = {
'dir_src'   : '%s/out' % os.environ['PLA1'],
'dir_dst'   : '%s/figs' % os.environ['PLA1'],
'id'        : (1, 2, 3, 4, 5, 6, 7, 8),
'Nm_slab'   : (64, 64, 64, 64, 256, 64, 64, 64),
'Nm_2d'     : (64, 128, 64, 64, 256, 128, 64, 64),
'eps_o'     : (3.33e-6, 1e-5, 3.33e-5, 4.64e-5, 1e-4, 3.33e-4, 4.64e-4, 1e-3),
'label'     : [], #('uno', 'dos', 'tres'),
}

for e in ps['eps_o']:
    ps['label'] += [ '$\epsilon: %1.2e$'%e ]

gp = sf.GralPlot(ps=ps, check=('eps_o',), check_all=True)
#gp = sf.GralPlot(ps=ps, check=None, check_all=False)
gp.do_checks()
gp.plot_kdiff()
gp.plot_errdy()
gp.plot_errEk(ylim=(1e-4, 10.))

#EOF
