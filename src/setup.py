#!/usr/bin/env ipython
from distutils.core import setup, Extension
from Cython.Build import cythonize

modname = 'cython_wrapper'

ext = Extension(
    name = modname,
    sources=[
        "%s.pyx" % modname, 
        "tt.cc", 
        "defs_turb.cc", 
        "funcs.cc",
        "general.cc", # para q sepa quien es 'scl'
        "odeintt.cc", 
        "stepperbs.cc",
    ],
    language="c++",
    #--- for debugging with 'gdb python'
    #extra_compile_args = ['-g'],
    #extra_link_args = ['-g'],
)

setup(
    name = modname,
    ext_modules = cythonize(ext)
)

#EOF
