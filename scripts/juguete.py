#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import os, sys
from os.path import isfile, isdir
import shared.funcs as sf

ps = {
'dir_src'   : '%s/out' % os.environ['PLA1'],
'dir_dst'   : '%s/figs' % os.environ['PLA1'],
'id'        : (9, 10, 11, 12, 5),
'Nm_slab'   : (64, 64, 128, 128, 256),
'Nm_2d'     : (64, 128, 128, 256, 256),
#'eps_o'     : (3.33e-6, 1e-5, 3.33e-5, 4.64e-5, 1e-4, 3.33e-4, 4.64e-4, 1e-3),
'label'     : [], #('uno', 'dos', 'tres'),
}

for Nms, Nm2d in zip(ps['Nm_slab'], ps['Nm_2d']):
    #ps['label'] += [ '$\epsilon: %1.2e$'%e ]
    ps['label'] += [ '$Nm^{slab}, Nm^{2d}: %d, %d$' % (Nms, Nm2d) ]

gp = sf.GralPlot(ps=ps, check=('eps_o',), check_all=True)
#gp = sf.GralPlot(ps=ps, check=None, check_all=False)
gp.do_checks()
gp.plot_kdiff()
gp.plot_errdy()
gp.plot_errEk(ylim=(1e-4, 10.))

#EOF
