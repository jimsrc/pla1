""" Script to smoke-test our Cython wrappers
"""
# Author: Gael Varoquaux
# License: BSD
import numpy as np

import cython_wrapper
a = cython_wrapper.py_compute(10)

print 'The array created is %s' % a
print 'It carries a reference to our deallocator: %s ' % a.base
np.testing.assert_allclose(a, np.arange(10))


b = cython_wrapper.py_compute_2d_ii(4, 4)
print " ---> 2d array: \n", b
print " ----> b.shape: ", b.shape
print 'It carries a reference to our deallocator: %s ' % b.base

raw_input()
