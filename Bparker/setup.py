#!/usr/bin/env ipython
from distutils.core import setup, Extension
from Cython.Build import cythonize

modname = 'Bparker'

ext = Extension(
    name = modname,
    sources=[
        "%s.pyx" % modname, 
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
