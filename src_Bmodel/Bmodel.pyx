# distutils: language = c++
# Author: Jimmy J. Masias-Meza
from libc.stdlib cimport free, malloc, calloc
from cpython cimport PyObject, Py_INCREF#, PyMem_Malloc, PyMem_Free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
#from cython.Utility.MemoryView import PyMem_New, PyMeM_Del # why doesn't work??
from cython.operator cimport dereference as deref
from libc.math cimport sqrt, sin, cos, M_PI, pow

# agregamos la clase wrapper
include "array_wrapper.pyx"

#--------------------------------------------
cdef class Bmodel(object):
    cdef PARAMS_TURB            *pt
    #cdef PARAMS                 *par
    cdef MODEL_TURB             *mt
    #cdef Odeint[StepperBS[rhs]] *bsode
    #cdef rhs                    *d
    #cdef Output[StepperBS[rhs]] *outbs
    cdef Doub tmax # t maximo absouto
    pdict = {}

    def __cinit__(self):
        self.mt = new MODEL_TURB()
        if self.mt is NULL:
            raise MemoryError()
        self.pt = new PARAMS_TURB()
        if self.pt is NULL:
            raise MemoryError()
        #self.outbs = new Output[StepperBS[rhs]]()
        #self.par = <PARAMS*> PyMem_Malloc(sizeof(PARAMS))

    def set_Bmodel(self, pdict, nB=0):
        """ inputs:
        - pdict   : diccionario de parametros de turbulencia
        - nB      : nro de B-realization
        """
        for nm in pdict.keys():
            self.pdict[nm] = pdict[nm]
        self._build_pturb() # objeto PARAMS_TURB
        self._build_par(nB=nB) # objeto MODEL_TURB
        #self.outbs.set_Bmodel(self.par)

    def _build_pturb(s):
        """
        build PARAMS_TURB object 'self.pt' from
        dictionary input 'self.pdict'
        """
        pd = s.pdict
        # parametros fisicos
        s.pt.Nm_slab = pd['Nm_slab']
        s.pt.Nm_2d   = pd['Nm_2d']
        s.pt.lmin_s  = pd['lmin_s']
        s.pt.lmax_s  = pd['lmax_s']
        s.pt.lmin_2d = pd['lmin_2d']
        s.pt.lmax_2d = pd['lmax_2d']
        s.pt.Lc_slab = pd['Lc_slab']
        s.pt.Lc_2d   = pd['xi']*pd['Lc_slab']
        s.pt.sigma_Bo_ratio = pd['sigma_Bo_ratio']
        s.pt.percent_slab = pd['ratio_slab']
        s.pt.percent_2d   = 1.0-pd['ratio_slab']
        # semillas
        s.pt.sem.slab[0] = pd['sem_slab0']
        s.pt.sem.slab[1] = pd['sem_slab1']
        s.pt.sem.slab[2] = pd['sem_slab2']
        s.pt.sem.two[0]  = pd['sem_two0']
        s.pt.sem.two[1]  = pd['sem_two1']

    def _build_par(s, nB=0):
        """ objeto PARAMS (todo):
        - aloco memoria para dB, B, etc..
        - defino modelo B con 'self.pt'
        - build dB spectra
        - fix B realization
        """
        ndim = 3
        #s.par.B       = PyMem_New(double, 3)
        # use 'PyMem_Malloc' instead of 'malloc'
        # src: https://docs.python.org/2/c-api/memory.html
        s.mt.B       = <Doub*> calloc(ndim, sizeof(Doub))
        s.mt.dB      = <Doub*> calloc(ndim, sizeof(Doub))
        s.mt.dB_SLAB = <Doub*> calloc(ndim, sizeof(Doub))
        s.mt.dB_2D   = <Doub*> calloc(ndim, sizeof(Doub))

        #NOTE: al 's.pt' lo construi en 'self._build_pturb()'
        s.mt.p_turb  = s.pt[0] #le paso todos los parametros! (MAGIA!!??!?)
        s.mt.p_turb.build_spectra()
        s.mt.fix_B_realization(nB=nB)

    def __dealloc__(self):
        """
        NOTE: cython will ignore
              a 'del self.pdict' because
              will consider it a read-only
        """
        if self.pt is not NULL:
            del self.pt
        if self.mt is not NULL:
            free(self.mt.B)
            free(self.mt.dB)
            free(self.mt.dB_SLAB)
            free(self.mt.dB_2D)
            del self.mt
            """PyMem_Free(self.mt)"""

    def Bxyz(self, xyz):
        cdef double pos[3]
        pos[0]=xyz[0]; pos[1]=xyz[1]; pos[2]=xyz[2]
        self.mt.calc_B(&(pos[0]))
        B = np.zeros(3)
        B[0]=self.mt.B[0]; B[1]=self.mt.B[1]; B[2]=self.mt.B[2]
        return B
#--------------------------------------------
#cpdef return_B(np.ndarray[np.float32_t, ndim=1, mode='c'] xyz):
#    """
#    input: position in spherical coords --> r [AU], th [deg], ph [deg]
#           as a numpy array.
#    output: 
#           tuple of Br, Bth, Bph # [Gauss]
#    """
#    AUincm = 1.5e13
#    cdef:
#        double pos[3], B[3]
#
#    pos[0] = xyz[0]*AUincm  # [cm] r, helioradius
#    pos[1] = xyz[1]*M_PI/180.  # [rad] th, co-latitude
#    #pos[2] = xyz[2]*AUincm  # [] phi --> B tiene simetria azimutal
#    calc_B(pos, B)
#    print B[0], B[1], B[2] # Bx, By, Bz [Gauss]
#    return B[0], B[1], B[2]


#EOF
