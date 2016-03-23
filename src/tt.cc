#ifndef TT_CC
#define TT_CC
#include "tt.h"

jim::jim(){
    a = 77;
    h1 = 0.077;
    aa = 0.9994;
}

/*
int run_orbit(PARAMS_TURB pturb){
	Output<StepperBS<rhs> > outbs;    // DES-COMENTAR
	outbs.set_Bmodel(pturb);            // DES-COMENTAR
	rhs d;		// object representing system-of-equations

    Odeint<StepperBS<rhs> > bsode(ystart,x1,x2,atol,rtol,h1,hmin,outbs,d,pturb,w_rank);	// inicializo el integrador para c/pla
    bsode.integrate();
}
*/


int run(double x1, double rigidity){
    double aux=2.0;
    scl.build(rigidity);
    printf(" scl.vel: %g\n", scl.vel);

	const Int nvar=6;		// nmbr of y-variables
	VecDoub ystart(nvar);		// allocate initial x, y[0, 1] values
    ystart[0] = 9.98;
    printf(" ys: %g\n", ystart[0]);


    return int(x1*aux);
}


int run3(double x1, double rigidity, PARAMS_TURB pt){
    double aux=2.0;
    printf(" @run3 ---> pt.nmod: %d\n", pt.n_modos);
    printf(" @run3 ---> pt.lambda_min: %g\n", pt.lambda_min);
    return int(x1*aux);
}

int run2(double x1, double rigidity, PARAMS_TURB* pt){
    double aux=2.0;
    return int(x1*aux);
}


#endif // TT_CC
