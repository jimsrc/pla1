#include "control.h"
#include "nr3.h"
#include "general.cc"

#include "defs_turb.h" // es incluido en funcs.h

#include "funcs.h"
#include "odeintt.h"	// "odeint.h"
#include "stepperbs.h"

extern ESCALAS scl;

int run_orbit(double *y_init, double x1, double x2){
	const Int nvar=6;		// nmbr of y-variables
	VecDoub ystart(nvar);		// allocate initial x, y[0, 1] values
    int i;
    // initialize init cond
    for(i=0; i<nvar; i++){
        ystart[i] = y_init[i];
    }
    double atol, rtol, h1, hmin;
    atol = rtol = 1e-5;
    h1=1e-10; hmin=0.0;

	//PARAMS par(fname_turb); 			// inputs para c/ region
    /*
	Output<StepperBS<rhs> > outbs;    // DES-COMENTAR
	outbs.set_Bmodel(par);            // DES-COMENTAR
	rhs d;		// object representing system-of-equations
    */

    //Odeint<StepperBS<rhs> > bsode(ystart,x1,x2,atol,rtol,h1,hmin,outbs,d,par,w_rank);	// inicializo el integrador para c/pla
    //bsode.integrate();
    
    return 0;
}

/*
int run_orbit(double *y, double x1, double x2){
    return 0;
}
*/
//EOF
