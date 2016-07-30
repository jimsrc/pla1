// defino mis estrucutras
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include "stdlib.h"

//-----------------CTES UNIVERSALES                                        
#define c               (3.0*1e10)              // [cm/s]                          
#define mo              (1.6726*1e-24)          // PROTON [gr]                     
#define E_reposo        938272013.0             // PROTON [eV]                     
#define q               (4.8032*1e-10)          // PROTON [statC]                  
#define AU_in_cm        1.5e13                  // [cm]
#define Vsw             (400*1e5)               // [cm/s]
#define omega           (2*M_PI/(27.27*86400))  // [1/s]
#define Bo              (5*1e-5)                // [G]  (1nT = 1e-5 G)
#define A               (Bo*AU_in_cm*AU_in_cm)  // [G*cm]

using namespace std;

typedef struct{ 	// pars_diff
        double Dperp;
        double a_perp;
        double b_perp;
        double Dparall;
        double a_parall;
        double b_parall;
} PARAMS_DIFF;

typedef struct{		//qcube_params
        int Ncubes;
        int Nr; 
        int Nth;
        int Nph;
        double dr; 
        double dth;
        double dph;
} PARAMS_QCUBE;

typedef struct{
	double rigidity;
	double gamma;
	double beta;
	double B[3];
	double Bmod;
} PARAMS_PLA;

void mostrar_parametros(double, double, double, double, double, PARAMS_DIFF *, double, int, string, double *, string); 
void actualiza(int, double *, double *); 
void setear_save_times(char *, string, double, double, int, double *, int &); 
void resetear(int, double *); 
void guarda_coords(int, double *, double *); 
void READ_INPUT(double &, string &, int &, double &, double &, PARAMS_DIFF *, double &, double &, double &, double &, string &, string);
double distrib(PARAMS_QCUBE *, double, double, double, double);
double ****alloc_quadri_cubo(PARAMS_QCUBE *); 
void acumula_en_distrib(PARAMS_QCUBE *, int, double *, double ****);
void OUTPUT(PARAMS_QCUBE *, double ****, ofstream&);
void integra_euler(long *SEM, double, PARAMS_PLA *pla, PARAMS_DIFF *par_diff, double *y, double dt, double *yout);
void particula_params(double RIGIDITY, PARAMS_PLA *pla);
void calc_B(double *, double, PARAMS_PLA *);
void set_seeds(int, int, long *, long *);
bool condic_contorno_interna(string, double, double, double &, double *, double *);
//
//---- FIN-------
