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
cdef extern from "exec.cc":
    int run_orbit(double* y_init, double x1, double x2)
    cdef cppclass jjj:
        jjj()
        int a

cdef extern from "tt.h":
    cpdef cppclass jim:
        jim()
        int a
        int b
        double aa

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

"""
cdef class jimmy:
    #wrapper around jimmy :-)
    cdef jim *_thisptr
    #cdef double aa 
    def __cinit__(self):
        self._thisptr = new jim()
        if self._thisptr == NULL:
            raise MemoryError()

    #def __init__(self):
    #   self.aa = self._thisptr.a

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    def set_aa(self, aa): # mejor usamos __set__()
        self._thisptr.aa = aa
        return self._thisptr.aa

    cpdef set_value(self, name, value):
        self._thisptr.aa = value

    property aa:
        def __get__(self):      return self._thisptr.aa
        def __set__(self, v):   self._thisptr.aa = v

    property var:
        def __get__(self):
            a = {
                'aa': self._thisptr.aa,
                'a': self._thisptr.a,
                'b': self._thisptr.b,
            }
            return a
"""

cdef class my_jjj:
    cdef jjj *_thisp
    def __cinit__(self):
        self._thisp = new jjj()
        if self._thisp == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisp!=NULL:
            del self._thisp

    property var:
        def __get__(self):
            v = {
                'n_modos'       : self._thisp.a,
            }
            return v

"""
cdef class PTurb:
    cdef PARAMS_TURB *_thisp
    def __cinit__(self):
        self._thisp = new PARAMS_TURB()
        if self._thisp == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisp!=NULL:
            del self._thisp

    property var:
        def __get__(self):
            a = {
                'n_modos'       : self.n_modos,
                'lambda_min'    : self.lambda_min,
            }
            return a
"""

cdef class ArrayWrapper:
    #We need to build an array-wrapper class to deallocate 
    #our array when the Python object is deleted.
    cdef void* data_ptr
    cdef int size

    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array

        This cannot be done in the constructor as it must recieve C-level
        arguments.

        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data            

        """
        self.data_ptr = data_ptr
        self.size = size

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_DOUBLE, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)



def py_run_orbit(np.ndarray[double, ndim=1, mode="c"] inp_array not None, double x1, double x2):
    """ 
    A wrapper around run_orbit() of exec.cc :-)

    inputs:
        inp_array   : input array
        x1          : start time
        x2          : end time
    """
    n = inp_array.shape[0]
    return run_orbit(&inp_array[0], x1, x2)

"""
def py_compute(int size):
    # Python binding of the 'compute' function in 'c_code.c' that does
    #    not copy the data allocated in C.
    #
    cdef double *array
    cdef np.ndarray ndarray
    # Call the C function
    array = compute(size)

    array_wrapper = ArrayWrapper()
    array_wrapper.set_data(size, <void*> array) 
    ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)

    return ndarray
"""
#EOF
