# distutils: language = c++

from libc.stdlib cimport free
from cpython.mem cimport PyMem_Free
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


#----------------------------------------------
cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size
    cdef bint survive

    cdef set_data(self, int size, void* data_ptr, survive=False):
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
        self.survive = survive

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
        references to the object are gone. 
        NOTE: In particular, for this project, I should not free
              the pointer since it must survive after 
              I "__get__()" it (that's the only motive why I use this 
              array-wrapper). So I must have 'self.survive=True'.
        """
        if not(self.survive) and (self.data_ptr is not NULL):
            print " ---> array_wrapper(1d): eliminando self.data_ptr @ ", self
            PyMem_Free(<void*>self.data_ptr)
            print " ---> eliminado ok!"



#----------------------------------------------
cdef class ArrayWrapper_2d:
    cdef void* data_ptr
    cdef int size, nx, ny
    cdef bint survive

    cdef set_data(self, int nx, int ny, void* data_ptr, survive=False):
        self.data_ptr = data_ptr
        self.nx = nx
        self.ny = ny
        self.survive = survive

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
        references to the object are gone. 
        NOTE: In particular, for this project, I should not free
              the pointer since it must survive after 
              I "__get__()" it (that's the only motive why I use this 
              array-wrapper). So I must have 'self.survive=True'.
        """
        if not(self.survive) and (self.data_ptr is not NULL):
            print " ---> array_wrapper(2d): eliminando self.data_ptr @ ", self
            PyMem_Free(<void*>self.data_ptr)
            print " ---> eliminado ok!"
#EOF
