""" Small Cython file to demonstrate the use of PyArray_SimpleNewFromData
in Cython to create an array from already allocated memory.

Cython enables mixing C-level calls and Python-level calls in the same
file with a Python-like syntax and easy type cohersion. See 
http://cython.org for more information
"""
# distutils: language = c++
# Author: Gael Varoquaux
# License: BSD

# Declare the prototype of the C function we are interested in calling
cdef extern from "c_code.cpp":
    double *compute(int size)
    double **compute_2d(int nx, int ny)
    double *compute_2d_ii(int nx, int ny)

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# We need to build an array-wrapper class to deallocate our array when
# the Python object is deleted.

cdef class ArrayWrapper:
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


cdef class ArrayWrapper_2d:
    cdef void* data_ptr
    cdef int size, nx, ny

    cdef set_data(self, int nx, int ny, void* data_ptr):
        self.data_ptr = data_ptr
        self.nx = nx
        self.ny = ny

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[2]
        shape[0] = <np.npy_intp> self.nx
        shape[1] = <np.npy_intp> self.ny
        # Create a 2D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(
                        nd = 2, 
                        dims = shape,
                        typenum = np.NPY_DOUBLE, 
                        data = self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)



def py_compute(int size):
    """ Python binding of the 'compute' function in 'c_code.c' that does
        not copy the data allocated in C.
    """
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


def py_compute_2d(int nx, int ny):
    cdef double **array
    cdef np.ndarray ndarray
    # Call the C function
    array = compute_2d(nx, ny)

    array_wrapper = ArrayWrapper_2d()
    array_wrapper.set_data(nx, ny, <void*> array) 
    ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)
    #Py_INCREF(array_wrapper)

    return ndarray


def py_compute_2d_ii(int nx, int ny):
    cdef double *array
    cdef np.ndarray ndarray
    # Call the C function
    array = compute_2d_ii(nx, ny)

    array_wrapper = ArrayWrapper_2d()
    array_wrapper.set_data(nx, ny, <void*> array) 
    ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)
    #Py_INCREF(array_wrapper)

    return ndarray
