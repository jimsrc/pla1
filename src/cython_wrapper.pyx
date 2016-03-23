""" Small Cython file to demonstrate the use of PyArray_SimpleNewFromData
in Cython to create an array from already allocated memory.

Cython enables mixing C-level calls and Python-level calls in the same
file with a Python-like syntax and easy type cohersion. See 
http://cython.org for more information
"""
# distutils: language = c++
# Author: Jimmy J.
#from libcpp.string cimport string

# Declare the prototype of the C function we are interested in calling
"""
cdef extern from "c_code.cpp":
    double *compute(int size)"""
#cdef extern from "exec.cc":
#    int run_orbit(double* y_init, double x1, double x2)
#    cdef cppclass jjj:
#        jjj()
#        int a

cdef extern from "tt.h":
    cpdef cppclass jim:
        jim()
        int a
        int b
        double aa
    
    int run(double x1)

cdef extern from "defs_turb.h":
    cpdef cppclass PARAMS_SEM:
        PARAMS_SEM()
        double a
    
    cpdef cppclass PARAMS_TURB:
        PARAMS_TURB()
        int n_modos
        double lambda_min

cdef extern from "funcs.h":
    double calc_gamma(double v);

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

# Import the Python-level symbols of numpy
#import numpy as np

# Import the C-level symbols of numpy
#cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
#np.import_array()


def run_py(double x1):
    return run(x1)


def c_gamma(double v):
    return calc_gamma(v)


cdef class psem:
    #wrapper around jimmy :-)
    cdef PARAMS_SEM *_thisp
    #cdef double aa 
    def __cinit__(self):
        self._thisp = new PARAMS_SEM()
        if self._thisp == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisp != NULL:
            del self._thisp

    property var:
        def __get__(self):
            v = {
                'a': self._thisp.a,
            }
            return v


cdef class pturb:
    #wrapper around jimmy :-)
    cdef PARAMS_TURB *_thisp
    #cdef double aa 
    def __cinit__(self):
        self._thisp = new PARAMS_TURB()
        if self._thisp == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisp != NULL:
            del self._thisp

    property var:
        def __get__(self):
            v = {
            'n_modos': self._thisp.n_modos,
            'lambda_min': self._thisp.lambda_min,
            }
            return v


cdef class jimmy:
    #wrapper around jimmy :-)
    cdef jim *_thisp
    #cdef double aa 
    def __cinit__(self):
        self._thisp = new jim()
        if self._thisp == NULL:
            raise MemoryError()

    #def __init__(self):
    #   self.aa = self._thisp.a

    def __dealloc__(self):
        if self._thisp != NULL:
            del self._thisp

    def set_aa(self, aa): # mejor usamos __set__()
        self._thisp.aa = aa
        return self._thisp.aa

    cpdef set_value(self, name, value):
        self._thisp.aa = value

    property aa:
        def __get__(self):      return self._thisp.aa
        def __set__(self, v):   self._thisp.aa = v

    property var:
        def __get__(self):
            a = {
                'aa': self._thisp.aa,
                'a': self._thisp.a,
                'b': self._thisp.b,
            }
            return a


#cdef class my_jjj:
#    cdef jjj *_thisp
#    def __cinit__(self):
#        self._thisp = new jjj()
#        if self._thisp == NULL:
#            raise MemoryError()
#
#    def __dealloc__(self):
#        if self._thisp!=NULL:
#            del self._thisp
#
#    property var:
#        def __get__(self):
#            v = {
#                'n_modos'       : self._thisp.a,
#            }
#            return v





#def py_run_orbit(np.ndarray[double, ndim=1, mode="c"] inp_array not None, double x1, double x2):
#    """ 
#    A wrapper around run_orbit() of exec.cc :-)
#
#    inputs:
#        inp_array   : input array
#        x1          : start time
#        x2          : end time
#    """
#    n = inp_array.shape[0]
#    return run_orbit(&inp_array[0], x1, x2)

#EOF
