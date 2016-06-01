#ifndef ODEINT_H
#define ODEINT_H
#include "nr3.h"
#include "funcs.h"
//#include "stepperbs.h"
//extern class PARAMS;

//---
template<class Stepper>
struct Odeint {
	static const Int MAXSTP=(150*50000); //MAXSTP=27*50000;
	Doub EPS;
	Int nok;
	Int nbad;
	Int nvar;			// dimension del problema: y[0], y[1], ..., y[nvar-1].
	Doub x1, x2, hmin;
	bool dense;			// True if dense output requested by out.
	VecDoub y, yold, dydx;		// (*)
	VecDoub &ystart;
	Output<Stepper> &out;
	typename Stepper::Dtype &derivs; // get the type of derivs from the stepper
	Stepper s;
	Int nstp;
	Doub x,h;
	Odeint(VecDoub_IO &ystartt,const Doub xx1,const Doub xx2,
		const Doub atol,const Doub rtol,const Doub h1,
		const Doub hminn,Output<Stepper> &outt,
        typename Stepper::Dtype &derivss, 
        PARAMS, int);
	void integrate();
	PARAMS par;
	//------------------- scattering stuff
	void save_history(void);
    #ifdef MONIT_SCATTERING
	void check_scattering(void);
    #endif //MONIT_SCATTERING 
	double mu_old, mu_new, Bmod, vmod, dtau;
	// ------------------ otros
	int wrank;

    #ifdef KILL_HANDLER
    static Odeint<Stepper>* _thisptr;
    void abort_mission(int signum); // remove current *.owned files
    #endif //KILL_HANDLER
// (*): en el constructor, paso las direcciones de memoria de al Stepper 's' para
// les haga las modificaciones q quiera.
};

#endif // ODEINT_H
//EOF
