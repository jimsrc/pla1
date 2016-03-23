
from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("cython_wrapper",
                sources=[
                    "cython_wrapper.pyx", 
                    "tt.cc", 
                    "defs_turb.cc", 
                    "funcs.cc",
                    "general.cc", # para q sepa quien es 'scl'
                    ],
                language="c++")

setup(name="cython_wrapper",
      ext_modules=cythonize(ext))

#EOF
