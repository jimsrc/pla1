#ifndef CONTROL
#define CONTROL

#define CYTHON     1 // stuff for cython handling
//#define KILL_HANDLER 1 // removes current *.owned files
#define MONIT_SCATTERING 1 // monitorea scattering de plas
#define MONIT_STEP 1 // monitorea el step
//#define BETA_CHECK 1 // checkea q siempre beta<1
#define NVAR    (6)

//--------------------------------------CTES UNIVERSALES
#define SIM_MAXSTP      (150*50000) // max nuber of steps in Odeint
#define clight          (3.0*1e10)              // [cm/s]
#define AU_in_cm        (1.5e13)                // [cm]
#define nT_in_G         (1.0*1e-5)              // [1G=1e5nT]

//----- operations
#define NORM(x,y,z)  (sqrt(x*x+y*y+z*z))


//--- some members better be public, when cython is ON
#ifdef CYTHON
    #define PRIVATE_OR_PUBLIC public
#else
    #define PRIVATE_OR_PUBLIC private
#endif //CYTHON


#endif //CONTROL
