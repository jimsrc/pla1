//
#ifndef GENERAL_CC
#define GENERAL_CC

#include "general.h"
#include "control.h"
//using namespace std;

void ESCALAS::build(const double RIGIDITY){
	// valores por defecto
	mo	= (1.6726*1e-24);	// [gr] masa PROTON
	Ereposo	= 938272013.0;		// [eV] energia de reposo PROTON 
	Z	= +1.;			// [e] carga (positiva) del proton 
	q	= (4.8032*1e-10);	// [statC] carga PROTON 
	B	= 5e-5;			// [G] 5nT en Gauss

    gamma   = pow(pow(Z*RIGIDITY/Ereposo, 2) + 1 , 0.5);
    beta    = pow(1. - 1./(gamma*gamma) , 0.5);
    vel     = beta*clight;                                       // [cm/s]
    Bo      = B;                                            // [G]
    wc      = q*Bo / (gamma*mo*clight);                          // [s^-1]
    rl      = vel / wc;                                     // [cm]
}

ESCALAS scl;    // global instance

#endif //GENERAL_CC
