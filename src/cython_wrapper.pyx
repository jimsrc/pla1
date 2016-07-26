# distutils: language = c++
# Author: Jimmy J.
#from libcpp.string cimport string

# Declare the prototype of the C function we are interested in calling

from libc.stdlib cimport free, malloc, calloc
from cpython cimport PyObject, Py_INCREF#, PyMem_Malloc, PyMem_Free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
#from cython.Utility.MemoryView import PyMem_New, PyMeM_Del # why doesn't work??
from cython.operator cimport dereference as deref
from libc.math cimport sqrt, sin, cos

# agregamos la clase wrapper
include "array_wrapper.pyx"


""" why this doesn't work??
cdef init_out(Output[StepperBS[rhs]] *op):
    o = out()
    o._thisp = op
    return o
"""


cdef class mgr:
    cdef Output[StepperBS[rhs]] *outbs
    cdef PARAMS_TURB            *pt
    cdef PARAMS                 *par
    cdef Odeint[StepperBS[rhs]] *bsode
    cdef rhs                    *d
    cdef Doub tmax # t maximo absouto
    pdict = {}

    def __cinit__(self):
        self.outbs = new Output[StepperBS[rhs]]()
        if self.outbs is NULL:
            raise MemoryError()

        #self.par = <PARAMS*> PyMem_Malloc(sizeof(PARAMS))
        self.par = new PARAMS()
        if self.par is NULL:
            raise MemoryError()

    def SetSim(self, **kargs):
        """
        **kargs: tmax, h1, hmin, mu, ph.
        """
        cdef Doub atol, rtol
        # estos parametros deben ir inmediatamente a la 
        # documentacion
        #rigidity    = kargs['rigidity']
        tmax        = kargs['tmax']
        h1          = kargs['FracGyroperiod']
        hmin        = kargs['hmin']
        mu          = kargs['mu']
        ph          = kargs['ph']
        rtol        = kargs['rtol']
        atol        = kargs['atol']
        #--- cond inic.
        cdef VecDoub *yini
        # posiciones a cero
        yini = new VecDoub(6, 0.0)
        # veloc inicial
        yini[0][1] = sqrt(1.-mu*mu)*cos(ph) # [1] vx
        yini[0][3] = sqrt(1.-mu*mu)*sin(ph) # [1] vy
        yini[0][5] = mu                     # [1] vz

        self.d = new rhs()

        cdef Doub x1, x2
        x1, x2 = 0.0, tmax #1e4
        self.bsode = new Odeint[StepperBS[rhs]](
            yini[0],x1,x2,atol,rtol,h1,hmin, 
            self.outbs[0],self.d[0],self.par[0], 0
        )
        #scl.build(rigidity) #(1e7)

    def RunSim(self):
        """ integramos una pla """
        self.outbs.tic()
        self.bsode.integrate()
        self.outbs.toc()

    def set_Bmodel(self, pdict, nB=0):
        """ inputs:
        - pdict   : diccionario de parametros de turbulencia
        - nB      : nro de B-realization
        """
        for nm in pdict.keys():
            self.pdict[nm] = pdict[nm]
        self._build_pturb() # objeto PARAMS_TURB
        self._build_par(nB=nB) # objeto PARAMS
        self.outbs.set_Bmodel(self.par)

    def build(self, str_timescale, nsave, tmaxHistTau, nHist, nThColl, i, j, dir_out):
        self.outbs.build(str_timescale, nsave, tmaxHistTau, nHist, nThColl, i, j, dir_out)

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
        s.par.B       = <Doub*> calloc(ndim, sizeof(Doub))
        s.par.dB      = <Doub*> calloc(ndim, sizeof(Doub))
        s.par.dB_SLAB = <Doub*> calloc(ndim, sizeof(Doub))
        s.par.dB_2D   = <Doub*> calloc(ndim, sizeof(Doub))

        #NOTE: al 's.pt' lo construi en 'self._build_pturb()'
        s.par.p_turb  = s.pt[0] #le paso todos los parametros! (MAGIA!!??!?)
        s.par.p_turb.build_spectra()
        s.par.fix_B_realization(nB=nB)

    def _build_pturb(s):
        """ build PARAMS_TURB object 'self.pt' from
        dictionary input 'self.pdict' """
        s.pt = new PARAMS_TURB()
        if s.pt is NULL:
            raise MemoryError()

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
        s.pt.percent_2d   = 1.0-pd['ratio_slab'] #pd['percent_2d']
        #s.pt.Bo      = pd['Bo']
        # semillas
        s.pt.sem.slab[0] = pd['sem_slab0']
        s.pt.sem.slab[1] = pd['sem_slab1']
        s.pt.sem.slab[2] = pd['sem_slab2']
        s.pt.sem.two[0]  = pd['sem_two0']
        s.pt.sem.two[1]  = pd['sem_two1']

    def __dealloc__(self):
        """NOTE: cython will ignore
              a 'del self.pdict' because
              will consider it a read-only"""
        if self.outbs is not NULL:
            del self.outbs
        if self.pt is not NULL:
            del self.pt
        if self.par is not NULL:
            free(self.par.B)
            free(self.par.dB)
            free(self.par.dB_SLAB)
            free(self.par.dB_2D)
            del self.par
            """PyMem_Free(self.par)"""

    def save2file(self):
        self.outbs.save2file()
    
    def Bxyz(self, xyz):
        cdef double pos[3]
        pos[0] = xyz[0]
        pos[1] = xyz[1]
        pos[2] = xyz[2]
        self.par.calc_B(&(pos[0]))
        B = np.zeros(3)
        B[0] = self.par.B[0]
        B[1] = self.par.B[1]
        B[2] = self.par.B[2]
        return B

    def sems(self):
        cdef long s0   = self.pt.sem.slab[0]
        cdef long two0 = self.pt.sem.two[0]
        return s0, two0

    property tsave:
        def __get__(self):
            cdef double *ptr
            cdef np.ndarray ndarray
            #n = self.outbs.xsave.size()
            n = self.outbs.count          # nro of saved times
            ptr = &(self.outbs.xsave[0])
            arrw = ArrayWrapper()
            arrw.set_data(n, <void*> ptr, survive=True)
            ndarray = np.array(arrw, copy=False)
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            return ndarray

    property ysave:
        def __get__(self):
            cdef double *ptr
            cdef np.ndarray ndarray
            cdef int nx, ny
            nx = 6                        # 'nvar' en c++
            #ny = self.outbs.count        # nro of saved times
            ny = self.outbs.ysave.ncols()
            arrw = ArrayWrapper_2d()

            ptr = &(self.outbs.ysave[0][0])
            arrw.set_data(nx, ny, <void*> ptr, survive=True)

            ndarray = np.array(arrw, copy=False)
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            #free(self.ptr)
            return ndarray

    property xyz:
        def __get__(self):
            cdef int n
            n = self.outbs.count
            v = nans((n, 3))
            v[:,0] = self.ysave[0,:n]   # [1]
            v[:,1] = self.ysave[2,:n]   # [1]
            v[:,2] = self.ysave[4,:n]   # [1]
            return v

    property vel:
        def __get__(self):
            n = self.outbs.count
            v = nans((n, 4))
            v[:,0] = self.ysave[1,:n]
            v[:,1] = self.ysave[3,:n]
            v[:,2] = self.ysave[5,:n]
            sum2   = np.empty(n, dtype=np.float64)
            sum2   = v[:,0]*v[:,0]+v[:,1]*v[:,1]+v[:,2]*v[:,2]
            v[:,3] = np.power(sum2, 0.5) # modulo
            return v

    property err:
        """ error de la energia total  """
        def __get__(self):
            """ retorna los errores de la energia total
            """
            cdef Doub vx, vy, vz, v, gamma
            cdef Int i, n
            n = self.outbs.count
            err = nans(n)
            for i in range(n):
                vx  = self.outbs.ysave[1][i]    # [1]
                vy  = self.outbs.ysave[3][i]    # [1]
                vz  = self.outbs.ysave[5][i]    # [1]
                v   = sqrt(vx*vx + vy*vy + vz*vz)
                err[i] = v-1.0 # veloc-relative-error
                #gamma   = calc_gamma(v);
                #err[i]  = gamma/scl.gamma - 1.0

            return err

    property mu:
        """ mu: cosine of pitch angle """
        def __get__(self):
            cdef double *ptr
            cdef np.ndarray ndarray
            n = self.outbs.count          # nro of saved times
            ptr = &(self.outbs.mu[0])
            arrw = ArrayWrapper()
            arrw.set_data(n, <void*> ptr, survive=True)
            ndarray = np.array(arrw, copy=False)
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            return ndarray

    property scl:
        def __get__(self):
            v = {
                'wc'    : scl.wc,
                'vel'   : scl.vel,
                'rl'    : scl.rl,
                'Bo'    : scl.Bo,
                'beta'  : scl.beta,
                'gamma' : scl.gamma,
            }
            return v

    property step_save:
        def __get__(self):
            cdef Doub *ptr
            cdef np.ndarray ndarray
            cdef int nx, ny
            nx = self.outbs.step_save.nrows()
            ny = self.outbs.step_save.ncols()
            arrw = ArrayWrapper_2d()
            ptr = &(self.outbs.step_save[0][0])
            #print self.outbs.HistStep[0][1]
            arrw.set_data(nx, ny, <void*> ptr, survive=True)

            ndarray = np.array(arrw, copy=False)
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            #free(self.ptr)
            return ndarray

    property Tau:
        def __get__(self):
            cdef Doub *ptr
            cdef np.ndarray ndarray
            cdef int nx, ny
            nx = self.outbs.Tau.nrows()
            ny = self.outbs.Tau.ncols()
            arrw = ArrayWrapper_2d()
            ptr = &(self.outbs.Tau[0][0])
            arrw.set_data(nx, ny, <void*> ptr, survive=True)
            ndarray = np.array(arrw, copy=False)
            #--- lo reducimos al trozo no-trivial
            ndarray = ndarray[:self.outbs.nreb,:]
            #------------------------------------
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            #free(self.ptr)
            return ndarray

    property gc_r:
        def __get__(self):
            cdef Doub *ptr
            cdef np.ndarray ndarray
            cdef int nx, ny
            nx = self.outbs.gc.r_gc.nrows()
            ny = self.outbs.gc.r_gc.ncols()
            arrw = ArrayWrapper_2d()
            ptr = &(self.outbs.gc.r_gc[0][0])
            arrw.set_data(nx, ny, <void*> ptr, survive=True)
            ndarray = np.array(arrw, copy=False)
            #--- lo reducimos al trozo no-trivial
            ndarray = ndarray[:self.outbs.gc.n,:]
            #------------------------------------
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            #free(self.ptr)
            return ndarray

    property gc_t:
        def __get__(self):
            #--- ahora el tiempo
            cdef double *ptr
            cdef np.ndarray ndarray
            #n = self.outbs.gc.n
            nx = self.outbs.gc.r_gc.nrows()
            ptr  = &(self.outbs.gc.t[0])
            arrw = ArrayWrapper()
            arrw.set_data(nx, <void*> ptr, survive=True)
            ndarray  = np.array(arrw, copy=False)
            #--- lo reducimos al trozo no-trivial
            ndarray  = ndarray[:self.outbs.gc.n]
            #------------------------------------
            ndarray.base = <PyObject*> arrw
            Py_INCREF(arrw)
            return ndarray


def nans(sh):
    return np.nan*np.ones(sh)



#def c_gamma(double v):
#    return calc_gamma(v)


#cdef void calc_Rlarmor(Doub Ek, Doub Bo, Doub *Rl):
#cpdef double calc_Rlarmor(Doub Ek, Doub Bo):
cpdef double calc_Rlarmor(Doub rigidity, Doub Bo):
    """
    input:
    Ek      : [eV] kinetic energy
    rigi..  : [V] rigidity
    Bo      : [G] magnetic field in Gauss
    output:
    Rl  : [cm] larmor radii
    """
    cdef:
        double q = (4.8032*1e-10) # [statC] carga PROTON
        double mo = 1.6726e-24 # [gr] masa PROTON
        double c = 3e10            # [cm/s] light speed
        double AU_in_cm = 1.5e13     # [cm]
        double E_reposo=938272013.0  # [eV] PROTON
        double beta, gamma, omg, v

    #rigidity = sqrt(Ek*Ek + 2.*Ek*E_reposo);
    #------------------------CALCULO DE GAMMA Y BETA
    gamma = pow(pow(rigidity/E_reposo,2) + 1. , 0.5)
    beta = pow(1. - 1/(gamma*gamma) , 0.5)
    #------------------------------CALCULO CICLOTRON
    omg = q * Bo / (gamma * mo * c)     # [s^-1]
    #---------------------------CALCULO RADIO LARMOR
    v   = beta * c              # [cm/s]
    #Rl[0]  = (v / omg) /AU_in_cm  # [AU]
    return (v / omg) # [cm]


#EOF
