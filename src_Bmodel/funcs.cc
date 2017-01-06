// incluyo librerias
# include <math.h>
# include <iostream>
# include <fstream>
# include <iomanip>
# include "stdlib.h"

//-----------------CTES UNIVERSALES                                        
#define clight          (3.0*1e10)              // [cm/s]   
#define mo              (1.6726*1e-24)          // PROTON [gr]  
#define E_reposo        938272013.0             // PROTON [eV]  
#define q               (4.8032*1e-10)          // PROTON [statC]
#define AU_in_cm        1.5e13                  // [cm]
#define Vsw             (400*1e5)               // [cm/s]
#define omega           (2*M_PI/(27.27*86400))  // [1/s]
#define Bo              (5*1e-5)                // [G]  (1nT = 1e-5 G)
#define A               (Bo*AU_in_cm*AU_in_cm)  // [G*cm]
#define SGN_B           (-1.0)                  // [1] HMF polarity

using namespace std;


// guarda en esfericas
void guarda_coords(int nsave, double *y, double *ysave){
    for(int i=0; i<3; i++){
        ysave[3*(nsave+1) + i] = y[i];      // para nsave=-1, guarda la cond inic
    }
}


double heaviside(double x){
    if(x<0)
        return 0.0;
    else if(x>=0)
        return 1.0;
}


void calc_gamma(double *y, double& GAMMA){
    double r, th;
    r   = y[0];
    th  = y[1];
    GAMMA = r * omega * sin(th) / Vsw;
}


void calc_relativ(double RIGIDITY, double &beta, double &gamma){
        gamma = pow(pow(RIGIDITY/E_reposo,2) + 1 , 0.5);
        beta = pow(1 - 1/(gamma*gamma) , 0.5);
}


void calc_B(double *y, double *B){
    double r, th, GAMMA;
    r   = y[0];
    th  = y[1];

    calc_gamma(y, GAMMA);

    // 'B' en esfericas
    B[0] = SGN_B*(1. - 2*heaviside(th - M_PI/2.)) * A / (r*r);   // ^r    [G]
    B[1] = 0.;                              // ^th    [G]
    B[2] = SGN_B*(1. - 2*heaviside(th - M_PI/2.)) * (-GAMMA) * A / (r*r); // ^ph   [G]

    //pla->Bmod   = pow(pla->B[0]*pla->B[0] + pla->B[1]*pla->B[1] + pla->B[2]*pla->B[2], 0.5);       //   [G]
}


void corrige_phi(double *y){
        // solo modifica 'ph' si cae fuera del intervalo [0, 2.*M_PI]
        double ph; 
        ph = y[2];
        if(ph > M_PI){
            ph = fmodf(ph, 2.*M_PI);
        }   
        else if(ph<0.){
            ph = 2.*M_PI + fmodf(ph, 2.*M_PI);
        }   
        y[2] = ph; 
}


void corrige_theta(double *y){
        // solo modifica 'th' si cae fuera de [0, M_PI]
        double th; 
        th      = y [1];
        if(th > M_PI){
            th = M_PI - fmodf(th, M_PI);
        }   
        else if(th < 0.){
            th = -th;
        }   
        y[1]    = th; 
}


/*
(*) le puse 1.97 xq si graficas esa expresion, es negativo entre 2 y 1.97. Luego es positivo a partir de 1.97 hasta 0.
La idea es q la contribucion sea positiva.
*/
