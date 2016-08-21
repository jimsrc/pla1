# distutils: language = c++
# Author: Jimmy J. Masias-Meza
import numpy as np
cimport numpy as np

# Numpy must be initialized. When using numpy from C or
# Cython you must _always_ do that, or you will have segfaults.
np.import_array()
#--- some datatypes 
#ctypedef np.int_t       Int
ctypedef np.float32_t   Float32
ctypedef np.ndarray     NDarray
#--- para q compile templates especificos
ctypedef double Doub
ctypedef int Int



cdef extern from "defs_turb.h":
    cpdef cppclass PARAMS_SEM:
        PARAMS_SEM()
        long slab[3]
        long two[2]

    cpdef cppclass PARAMS_TURB:
        PARAMS_TURB()
        void build_spectra()
        Int Nm_slab, Nm_2d
        Doub lmin_s, lmax_s
        Doub lmin_2d, lmax_2d
        Doub Lc_slab, Lc_2d
        Doub sigma_Bo_ratio;
        Doub percent_slab, percent_2d;
        Doub gS, g2D # potencia espectral slab/2D
        Doub sigma_S    # intensidad slab
        Doub sigma_2D   # intensidad 2D
        PARAMS_SEM sem  # semillas

    cpdef cppclass MODEL_TURB:
        MODEL_TURB()
        void fix_B_realization(const int nB)
        void calc_B(const Doub *pos)
        Doub *B
        Doub *dB
        Doub *dB_2D
        Doub *dB_SLAB
        PARAMS_TURB p_turb

#EOF
