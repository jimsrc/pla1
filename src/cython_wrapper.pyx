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
        PARAMS_SEM s
    
    int run(double x1, double rigidity)
    int run2(double x1, double rigidity, PARAMS_TURB* ptt)
    int run3(double x1, double rigidity, PARAMS_TURB pt);

cdef extern from "defs_turb.h":
    cpdef cppclass PARAMS_SEM:
        PARAMS_SEM()
        double a
    
    cpdef cppclass PARAMS_TURB:
        PARAMS_TURB()
        int n_modos
        double lambda_min
        PARAMS_SEM sem

    cpdef cppclass MODEL_TURB:
        MODEL_TURB()
        double *B
        double *dB
        double *dB_2D
        double *dB_SLAB
        PARAMS_TURB p_turb
        void calc_B(const double *)

cdef extern from "funcs.h":
    cpdef cppclass PARAMS:
        PARAMS()
        double *B
        double *dB
        double *dB_2D
        double *dB_SLAB
        PARAMS_TURB p_turb
        void calc_B(const double *)
        void test()

    double calc_gamma(double v);

    cdef cppclass Output[T]:
        Output()
        int nvar
        #void set_Bmodel(PARAMS par)

cdef extern from "odeintt.h":
    cdef cppclass Odeint[T]:
        Odeint()
        int nok

"""cdef extern from "stepperbs.h":
    cdef struct StepperBase
    cdef struct StepperBS"""

"""
cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        vector() except +
        vector(vector&) except +
        vector(size_t) except +
        vector(size_t, T&) except +
        T& operator[](size_t)
        void clear()
        void push_back(T&)
        size_t size()
"""

from libc.stdlib cimport free, malloc, calloc
from cpython cimport PyObject, Py_INCREF
from cython.operator cimport dereference as deref

# Import the Python-level symbols of numpy
#import numpy as np

# Import the C-level symbols of numpy
#cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
#np.import_array()


def run_py(double x1, double rig):
    cdef PARAMS_TURB *pt
    pt = new PARAMS_TURB()
    pt.n_modos = 71+3
    pt.lambda_min = 0.337
    print " --> pt.lam: ", pt.lambda_min
    run2(x1, rig, pt)
    #run3(x1, rig, deref(pt)) # deref(pt)==pt[0]
    run3(x1, rig, pt[0]) # deref(pt)==pt[0]

    #--- creamos un array :)
    cdef double* array;
    array = <double*> calloc(5, sizeof(double)) 
    print " array... ", array[0]

    #--- puntero al modelo
    cdef MODEL_TURB *mt
    mt = new MODEL_TURB()
    # alocamos los campos
    mt.B        = <double*> calloc(3, sizeof(double))
    mt.dB       = <double*> calloc(3, sizeof(double))
    mt.dB_SLAB  = <double*> calloc(3, sizeof(double))
    mt.dB_2D    = <double*> calloc(3, sizeof(double))
    # lo sgte FUNCIONA OK! :O :O ......... :O!!!!
    mt.p_turb   = pt[0] # que SUCEDE AQUI???
    print " --->>> mt.pturb: ", mt.p_turb.n_modos
    print " mt.B[0]: ", mt.B[2]

    #--- puntero a PARAMS (todo)
    cdef PARAMS *par
    par = new PARAMS()
    par.B        = <double*> calloc(3, sizeof(double))
    par.dB       = <double*> calloc(3, sizeof(double))
    par.dB_SLAB  = <double*> calloc(3, sizeof(double))
    par.dB_2D    = <double*> calloc(3, sizeof(double))
    par.p_turb   = pt[0] # que SUCEDE AQUI???
    print " --->>> mt.pturb: ", par.p_turb.n_modos
    par.test()

    #cdef vector[int] v
    #print "v: ", v.size()
    #cdef Odeint *ode
    #ode = new Odeint()
    #ode[0].nok = 666

    #--- free memory
    free(par.B);         free(mt.B)
    free(par.dB);        free(mt.dB)
    free(par.dB_SLAB);   free(mt.dB_SLAB)
    free(par.dB_2D);     free(mt.dB_2D)
    free(par);           free(mt)
    free(pt)
    # return
    return run(x1, rig)


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
