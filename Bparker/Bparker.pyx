# distutils: language = c++
# Author: Jimmy J.

from libc.stdlib cimport free, malloc, calloc
from cpython cimport PyObject, Py_INCREF#, PyMem_Malloc, PyMem_Free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
#from cython.Utility.MemoryView import PyMem_New, PyMeM_Del # why doesn't work??
from cython.operator cimport dereference as deref
from libc.math cimport sqrt, sin, cos, M_PI
import numpy as np
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()
#--- some datatypes 
ctypedef np.int_t       Int
ctypedef np.float32_t   Float32
ctypedef np.ndarray     NDarray

# agregamos la clase wrapper
include "array_wrapper.pyx"


cdef extern from "funcs.cc":
    void calc_B(double *y, double *B);


cpdef return_B(np.ndarray[np.float32_t, ndim=1, mode='c'] xyz):
    """
    input: position in spherical coords --> r [AU], th [deg], ph [deg]
           as a numpy array.
    output: 
           tuple of Br, Bth, Bph # [Gauss]
    """
    AUincm = 1.5e13
    cdef:
        double pos[3], B[3]

    pos[0] = xyz[0]*AUincm  # [cm] r, helioradius
    pos[1] = xyz[1]*M_PI/180.  # [rad] th, co-latitude
    #pos[2] = xyz[2]*AUincm  # [] phi --> B tiene simetria azimutal
    calc_B(pos, B)
    print B[0], B[1], B[2] # Bx, By, Bz [Gauss]
    return B[0], B[1], B[2]

