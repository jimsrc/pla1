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
    pdict = {}

    def __cinit__(self):
        self.outbs = new Output[StepperBS[rhs]]()
        if self.outbs is NULL:
            raise MemoryError()


    def runsim(self, **kargs):
        """
        **kargs:
            rigidity:
            tmax
            h1 
            hmin
            mu
            ph
        """
        cdef Doub atol, rtol
        # estos parametros deben ir inmediatamente a la 
        # documentacion
        rigidity    = kargs['rigidity']
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
        #atol = rtol = 1e-5

        cdef rhs *d
        d = new rhs()

        #cdef Odeint[StepperBS[rhs]] *bsode
        cdef Doub x1, x2
        x1, x2 = 0.0, tmax #1e4
        self.bsode = new Odeint[StepperBS[rhs]](
                yini[0],x1,x2,atol,rtol,h1,hmin, 
                self.outbs[0],d[0],self.par[0], 0
                )
        scl.build(rigidity) #(1e7)

        #--- integramos una pla
        self.outbs.tic()
        self.bsode.integrate()
        self.outbs.toc()

        print " ==> nrows ", self.outbs.ysave.nrows()
        print " ==> ncols ", self.outbs.ysave.ncols()
    

    def set_sim(self, **kargs):
        """
        **kargs:
            rigidity:
            tmax
            h1 
            hmin
            mu
            ph
        """
        cdef Doub atol, rtol
        # estos parametros deben ir inmediatamente a la 
        # documentacion
        rigidity    = kargs['rigidity']
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
        #atol = rtol = 1e-5

        cdef rhs *d
        d = new rhs()

        #cdef Odeint[StepperBS[rhs]] *bsode
        cdef Doub x1, x2
        x1, x2 = 0.0, tmax #1e4
        self.bsode = new Odeint[StepperBS[rhs]](
                yini[0],x1,x2,atol,rtol,h1,hmin, 
                self.outbs[0],d[0],self.par[0], 0
                )
        scl.build(rigidity) #(1e7)
        """
        #--- integramos una pla
        self.outbs.tic()
        self.bsode.integrate()
        self.outbs.toc()"""


    def runsim_partial(self, x2=None):
        if x2 is None: # chekear q sea tipo 'Doub'
            print " x2 must have a value!"
            raise SystemExit

        self.bsode.x2 = x2  # <double> solo modificamos el 'x2' del 'Odeint()'
        #--- integramos una pla
        self.outbs.tic()
        # chekear q el ystart sea el ultimo de la ultima "corrida"
        self.bsode.integrate()
        self.outbs.toc()
        """
        template<class Stepper>
        void Odeint<Stepper>::integrate(Doub xx2) {
            x2 = xx2;
            integrate();
        }"""


    def set_Bmodel(self, pdict, nB=0):
        """ inputs:
        - pdict   : diccionario de parametros de turbulencia
        - nB      : nro de B-realization
        """
        for nm in pdict.keys():
            self.pdict[nm] = pdict[nm]
        #self.pdict = pdict # no esta permitido cambiar el puntero!
        self.build_pturb() # objeto PARAMS_TURB
        self.build_par(nB=nB) # objeto PARAMS
        self.outbs.set_Bmodel(self.par[0])


    def build(self, str_timescale, nsave, tmaxHistTau, nHist, i, j, dir_out):
        self.outbs.build(str_timescale, nsave, tmaxHistTau, nHist, i, j, dir_out)


    def build_par(s, nB=0):
        """ objeto PARAMS (todo):
        - aloco memoria para dB, B, etc..
        - defino modelo B con 'self.pt'
        - build dB specta
        - fix B realization
        """
        #s.par = <PARAMS*> PyMem_Malloc(sizeof(PARAMS))
        s.par = new PARAMS()
        if s.par is NULL:
            raise MemoryError()

        ndim = 3
        #s.par.B       = PyMem_New(double, 3)
        # use 'PyMem_Malloc' instead of 'malloc'
        # src: https://docs.python.org/2/c-api/memory.html
        s.par.B       = <double*> calloc(ndim, sizeof(double))
        s.par.dB      = <double*> calloc(ndim, sizeof(double))
        s.par.dB_SLAB = <double*> calloc(ndim, sizeof(double))
        s.par.dB_2D   = <double*> calloc(ndim, sizeof(double))
        
        s.par.p_turb  = s.pt[0] #le paso todos los parametros! (MAGIA!!??!?)
        s.par.p_turb.build_spectra()
        s.par.fix_B_realization(nB=nB)


    def build_pturb(s):
        """ build PARAMS_TURB object 'self.pt' from
        dictionary input 'self.pdict' """
        s.pt = new PARAMS_TURB()
        if s.pt is NULL:
            raise MemoryError()

        pd = s.pdict
        # parametros fisicos
        s.pt.n_modos = pd['n_modos']
        s.pt.lambda_min = pd['lambda_min']
        s.pt.lambda_max = pd['lambda_max']
        s.pt.Lc_slab = pd['Lc_slab']
        s.pt.Lc_2d = pd['Lc_2d']
        s.pt.sigma_Bo_ratio = pd['sigma_Bo_ratio']
        s.pt.percent_slab = pd['percent_slab']
        s.pt.percent_2d = pd['percent_2d']
        s.pt.Bo = pd['Bo']
        # semillas
        s.pt.sem.slab[0] = pd['sem_slab0']
        s.pt.sem.slab[1] = pd['sem_slab1']
        s.pt.sem.slab[2] = pd['sem_slab2']
        s.pt.sem.two[0]  = pd['sem_two0']
        s.pt.sem.two[1]  = pd['sem_two1']


    """ for some reason, this compiles, but doesn't 
        run properly :(
    def set_all(self, pdict, nB, pother):
        if None in (pdict, nB, pother):
            print ' ---> argumento None en set_all(..)!!'
            raise SystemExit

        self.set_Bmodel(pdict=pdict, nB=nB)
        self.build(**pother)
    """


    def __dealloc__(self):
        #NOTE: cython will ignore
        #      a 'del self.pdict' because
        #      will consider it a read-only
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
            #PyMem_Free(self.par)


    def save2file(self):
        self.outbs.save2file()
    

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
                v 	= sqrt(vx*vx + vy*vy + vz*vz)
                gamma	= calc_gamma(v);
                err[i]	= gamma/scl.gamma - 1.0

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



def nans(sh):
    return np.nan*np.ones(sh)



def c_gamma(double v):
    return calc_gamma(v)


#EOF
