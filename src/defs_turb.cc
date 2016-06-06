#ifndef DEFS_TURB_CC
#define DEFS_TURB_CC


#include "control.h"
#include "general.h"
#include "defs_turb.h"
#include "funcs.h"


// declaration of definition in general.cc
//extern ESCALAS scl;


//=========================================================================

/*---------------------DEFs para "rann0" ---------------*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float rann0(long &idum){        // lo nombre diferente a "ran0" para evitar 
  long k;                       // posibles conflictos con las librerias
  float ans;
  idum ^= MASK;                   // XORing with MASK allows use of 0 and
  k = idum/IQ;                    //     other simple bit patterns for idum.
  idum = IA * (idum-k*IQ) - IR*k; // Compute idum = (IA*idum) % IM without
  if (idum < 0) idum += IM;       //     overflows by Schrage's method.
  ans = AM * idum;                // Convert idum to a floating result.
  idum ^= MASK;                   // Unmask before return.
  return ans;
}
/*-------------------------------------------------------*/
/*FASES::~FASES(){
    cout << " xxx destruyendo FASES." << endl;
}*/

void FASES::build(Int nm_slab, Int nm_2d, PARAMS_SEM *sem){ // constructor
    cout << " ...construyendo FASES." << endl;
    Nm_slab = nm_slab;
    Nm_2d   = nm_2d;
    //n_modos = n;
    phi_s   = new double[Nm_slab];
    a_s     = new double[Nm_slab];
    b_s     = new double[Nm_slab];
    phi_2d  = new double[Nm_2d];
    b_2d    = new double[Nm_2d];

    construir_fases_random(*sem); // le mando *sem para q no me cambie los valores originales
}

void FASES::construir_fases_random(PARAMS_SEM sem){
    //------- preparamos las semillas
    // "Aleatorizar" las semillas.
    // Esto es xq no queremos el 1er valor
    // random generado, pues x alguna razon, son 
    // proporcionales a las semillas, en el rango 
    // de semillas:(0, 2e+9).
    for(int i=0;i<3;i++) rann0(sem.slab[i]);
    for(int i=0;i<2;i++) rann0(sem.two[i]);
    //------------------------------
    
    for(int i=0; i<Nm_slab; i++){
        // fases para SLAB
        phi_s[i]    = 2.*M_PI* rann0(sem.slab[0]);
        a_s[i]      = 2.*M_PI* rann0(sem.slab[1]);
        b_s[i]      = 2.*M_PI* rann0(sem.slab[2]);
    }

    for(int i=0; i<Nm_2d; i++){
        //cout << " ----> semtwo1: "<< sem.two[1] << endl;
        // fases para 2D
        phi_2d[i]   = 2.*M_PI* rann0(sem.two[0]);
        b_2d[i]     = 2.*M_PI* rann0(sem.two[1]);
        //printf(" b2d(%d): %f\n", i, b_2d[i]);
    }
}




//---------------------------------------- class PARAMS_TURB
//PARAMS_TURB::PARAMS_TURB(): n_modos(10), lambda_min(4.567) { // esto funciona c/cython!!!!
PARAMS_TURB::PARAMS_TURB(){
    cout << " ...construyendo PARAMS_TURB *sin* input." << endl;
}
PARAMS_TURB::PARAMS_TURB(string fname_input){
    cout << " ...construyendo PARAMS_TURB *con* input: " << fname_input << endl << endl;
    build(fname_input);
}


/*PARAMS_TURB::~PARAMS_TURB(){
    cout << " xxx destruyendo PARAMS_TURB: " << FNAME_INPUT << endl;
}*/


void PARAMS_TURB::report(void){
    cerr << " REPORTE DE OBJETO PARAMS_TURB:-----------------------------" << endl;
    cerr << " Nm_slab:      " << Nm_slab        << endl;
    cerr << " Nm_2d:      " << Nm_2d          << endl;
    cerr << " lmin_s [AU]:  " << lmin_s / AU_in_cm  << endl;
    cerr << " lmax_s [AU]:  " << lmax_s / AU_in_cm  << endl;
    cerr << " lmin_2d [AU]:  " << lmin_2d / AU_in_cm  << endl;
    cerr << " lmax_2d [AU]:  " << lmax_2d / AU_in_cm  << endl;
    cerr << " Lc_slab [AU]:     " << Lc_slab / AU_in_cm         << endl;
    cerr << " Lc_2d   [AU]:     " << Lc_2d / AU_in_cm       << endl;
    cerr << " Bo [nT]:      " << Bo / nT_in_G       << endl;
    cerr << " sigma_Bo_ratio [1]:   " << sigma_Bo_ratio         << endl;
    cerr << " percent_slab [%]: " << 100.*percent_slab      << endl;
    cerr << " percent_2d [%]:   " << 100.*percent_2d        << endl;
    cerr << " gS [1]: " << gS << endl;
    cerr << " g2D [1]: " << g2D << endl;
    cerr << " -----------------------------------------------------------" << endl;
}


void PARAMS_TURB::build(string fname_input){
    FNAME_INPUT = fname_input;
    read_params(fname_input); // setea Nm_slab, Lc_slab, Lc_2d, lambda_min, etc...
    build_spectra();
}


/* Deben estar definidos:
 * Nm_slab          [1]
 * Nm_2d            [1]
 * sigma_Bo_ratio   [1]
 * Bo               [G]
 * percent_slab     [1] fraction
 * percent_2d       [1] fraction
 * lmin_s,  lmax_s  [cm] escalas de turb slab
 * lmin_2d, lmax_2d [cm] escalas de turb 2d
 * Lc_slab          [cm]
 * Lc_2d            [cm]
 * sem.slab[]       [1] tres semillas
 * sem.two[]        [1] dos semillas
*/
void PARAMS_TURB::build_spectra(){
    dk_s    = new Doub[Nm_slab];
    dk_2d   = new Doub[Nm_2d];
    k_s     = new Doub[Nm_slab];
    k_2d    = new Doub[Nm_2d];
    Bk_SLAB = new Doub[Nm_slab];
    Bk_2D   = new Doub[Nm_2d];

    fases.build(Nm_slab, Nm_2d, &sem);
    build_sigmas();
    build_k_and_dk();
    build_Bk_SLAB();
    build_Bk_2D();
    report();
}


void PARAMS_TURB::read_params(string fname_input){
    string dummy;
    ifstream filein(fname_input.c_str());
    if (!filein.good()) {
        cout << " problema al abrir " << fname_input << endl;
        exit(1);
    }

    // parametros del modelo
    filein >> Nm_slab     >> dummy; // [1] (entero) nro modos Slab
    filein >> Nm_2d       >> dummy; // [1] (entero) nro modos 2D
    filein >> lmin_s  >> dummy; // [AU] escala max Slab
    filein >> lmax_s  >> dummy; // [AU] escala min Slab
    filein >> lmin_2d >> dummy; // [AU] escala min 2D
    filein >> lmax_2d >> dummy; // [AU] escala max 2D
    filein >> Lc_slab     >> dummy; // [AU] longitud de correlacion Lc, SLAB 
    filein >> Lc_2d       >> dummy; // [AU] longitud de correlacion Lc, 2D
    filein >> Bo          >> dummy; // [G] campo medio Bo
    filein >> sigma_Bo_ratio >> dummy; // [1] sigma^2 = (sigma_Bo_ratio) * Bo^2 
    filein >> percent_slab   >> dummy; // [fraccion] sigma_SLAB^2 = percent_slab * sigma^2
    filein >> percent_2d     >> dummy; // [fraccion] sigma_2D^2 = percent_2D * sigma^2a

    // semillas
    filein >> sem.slab[0]       >> dummy;
    filein >> sem.slab[1]       >> dummy;
    filein >> sem.slab[2]       >> dummy;
    filein >> sem.two[0]        >> dummy;
    filein >> sem.two[1]        >> dummy;
    // correcc a unidades fisicas
    Lc_slab     *= AU_in_cm;        // [cm]
    Lc_2d       *= AU_in_cm;        // [cm]
    lmax_s      *= AU_in_cm;        // [cm]
    lmin_s      *= AU_in_cm;        // [cm]
    lmax_2d     *= AU_in_cm;        // [cm]
    lmin_2d     *= AU_in_cm;        // [cm]

}


void PARAMS_TURB::build_sigmas(){
    Doub sigma;
    sigma    = sqrt(sigma_Bo_ratio) * Bo;  // [G]
    sigma_S  = sqrt(percent_slab) * sigma; // [G]
    sigma_2D = sqrt(percent_2d) * sigma;   // [G]
}


void PARAMS_TURB::build_k_and_dk(){
    Doub kmin, kmax;

    //--- wavevectors slab
    kmin = (2. * M_PI) / lmax_s;     // [cm^-1]
    kmax = (2. * M_PI) / lmin_s;     // [cm^-1]
    for(int i=0; i<Nm_slab; i++){
        k_s[i]  = kmin * pow(kmax/kmin, 1.*i/(Nm_slab-1.));       // [cm^-1]
        dk_s[i] = k_s[i] * (pow(kmax/kmin, 1./(Nm_slab-1)) - 1.); // [cm^-1]
    }

    //--- wavevectors 2d
    kmin = (2. * M_PI) / lmax_2d;     // [cm^-1]
    kmax = (2. * M_PI) / lmin_2d;     // [cm^-1]
    for(int i=0; i<Nm_2d; i++){
        k_2d[i]  = kmin * pow(kmax/kmin, 1.*i/(Nm_2d-1.));       // [cm^-1]
        dk_2d[i] = k_2d[i] * (pow(kmax/kmin, 1./(Nm_2d-1)) - 1.);   // [cm^-1]
    }
}


void PARAMS_TURB::build_Bk_SLAB(){
    Int i;
    Doub DENOMINADOR=0.0, FACTOR;
    for(i=0; i<Nm_slab; i++){
        DENOMINADOR += dk_s[i] / (1. + pow(k_s[i]*Lc_slab, gS));
    }
    for(i=0; i<Nm_slab; i++){
        FACTOR = dk_s[i] / (1. + pow(k_s[i]*Lc_slab, gS)) / DENOMINADOR;
        Bk_SLAB[i] = sigma_S * pow(FACTOR, 0.5);            // [G]
    }
}


void PARAMS_TURB::build_Bk_2D(){
    Int i;
    Doub DENOMINADOR=0.0, FACTOR, dV;
    for(i=0; i<Nm_2d; i++){
        dV = 2.*M_PI*k_2d[i]*dk_2d[i];
        DENOMINADOR += dV / (1. + pow(k_2d[i]*Lc_2d, g2D));
    }
    for(i=0; i<Nm_2d; i++){
        dV = 2.*M_PI*k_2d[i]*dk_2d[i];
        FACTOR = dV / (1. + pow(k_2d[i]*Lc_2d, g2D)) / DENOMINADOR;
        Bk_2D[i] = sigma_2D * pow(FACTOR, 0.5);             // [G]
    }
}




//-------------------------------------- class MODEL_TURB

MODEL_TURB::MODEL_TURB(string fname_input) {  // constructor
    build(fname_input);
}   


// TODO: hacer a dB, dB_SLAB, dB_2D arrays con direcciones de memo 
//       continuas! Usar malloc() en vez de 'new double'. Propuesta:
//       fields  = double[3*4];  // esto en "defs_turb.h"
//       // esto en esta rutina
//       dB      = fields
//       dB_SLAB = fields + 3
//       dB_2D   = fields + 6
//       B       = fields + 9
//
void MODEL_TURB::build(string fname_input){ // construyo leyendo los params desde archivo
    FNAME_INPUT = fname_input;
    cout <<endl<< " ...construyendo MODEL_TURB: " << FNAME_INPUT << endl<<endl;

    dB      = new Doub[3];
    dB_SLAB = new Doub[3];
    dB_2D   = new Doub[3];
    B       = new Doub[3];

    // NOTA: en particular, debemos asegurarnos de q la
    // componente z del campo turb sea nula (y asi se quedara
    // en el resto del codigo):
    for(int i=0; i<3; i++){
        dB[i]       = 0.0;
        dB_SLAB[i]  = 0.0;
        dB_2D[i]    = 0.0;
    }

    // seteo: semillas, los modos, las fases, el espectro Fourier, etc.
    p_turb.build(fname_input);          // inicializo parametros/semillas del modelo
}


void MODEL_TURB::fix_B_realization(const int nBrz){
    long seed = 100*nBrz;   // semilla queda en funcion del nro de realizacion
    rann0(seed);rann0(seed);//para "aleatorizar" a 'seed'
    rann0(seed); p_turb.sem.slab[0] = seed;
    rann0(seed); p_turb.sem.slab[1] = seed;
    rann0(seed); p_turb.sem.slab[2] = seed;
    rann0(seed); p_turb.sem.two[0]  = seed;
    rann0(seed); p_turb.sem.two[1]  = seed;

    p_turb.fases.construir_fases_random(p_turb.sem);
}


// TODO: convertir estas variables COSCOS, phi, etc, en miembros 
//       privados de MODEL_TURB
// TODO: tal vez usar sin(x)=sqrt(1-cos(x)^2) resulte mas barato 
//       computacionalmente? (cuidado q sqrt() SIEMPRE DA POSITIVO)
// TODO: para evitar tener a 'p_turb' como publico, podriamos poner un
//       macro para q sea publico solo con Cython.
// TODO: servira de algo hacer copias de p_turb.fases.a_s[], etc para
//       accesar mas rapido a esos valores aqui?
//       Dichas copias las haria en MODEL_TURB::build(string)
void MODEL_TURB::calc_dB_SLAB(const Doub *pos){
    Doub COSCOS, SINSIN, FACTOR_X, FACTOR_Y, k;
    Doub b, a, phi;                     // fases random Slab
    dB_SLAB[0] = 0.; dB_SLAB[1]=0.;     // reseteo el vector en ceros

    for(Int i=0; i<p_turb.Nm_slab; i++){
        a           = p_turb.fases.a_s[i];
        b           = p_turb.fases.b_s[i];
        phi         = p_turb.fases.phi_s[i];
        k           = p_turb.k_s[i];
        COSCOS      = cos(k*pos[2] + b) * cos(a);
        SINSIN      = sin(k*pos[2] + b) * sin(a);
        FACTOR_X    = COSCOS * cos(phi) + SINSIN * sin(phi);
        FACTOR_Y    = COSCOS * sin(phi) - SINSIN * cos(phi);

        dB_SLAB[0] += p_turb.Bk_SLAB[i] * FACTOR_X;
        dB_SLAB[1] += p_turb.Bk_SLAB[i] * FACTOR_Y;
    }
}


// TODO: convertir declaraciones en variables privadas!, y 
//       probar "const Doub const *pos".
// TODO: esta calculando 'sin(FACTOR)' DOS VECES!!!
//       Borrar 'FACTOR', y usemos:
//       Doub sinFACTOR;
//       sinFACTOR = sin( k*(pos[0]...) );
//       dB_2D[0] += ... * sinFACTOR;
void MODEL_TURB::calc_dB_2D(const Doub *pos){
    Doub FACTOR, k;
    Doub phi, b;                 // fases random 2D
    dB_2D[0]=0.0; dB_2D[1]=0.0;  // reseteo el vector en ceros

    for(Int i=0; i<p_turb.Nm_2d; i++){
        b        = p_turb.fases.b_2d[i];
        phi      = p_turb.fases.phi_2d[i];
        k        = p_turb.k_2d[i];
        FACTOR   = k*(pos[0]*cos(phi) + pos[1]*sin(phi)) + b;
        dB_2D[0] += p_turb.Bk_2D[i] * sin(phi) * sin(FACTOR);
        dB_2D[1] +=-p_turb.Bk_2D[i] * cos(phi) * sin(FACTOR);
    }
}


// TODO: estos if(..>0.0) se pueden evitar con una variable const, asi:
//       --------------
//       const void calc_dB_some;
//       const bool a = check_pure_model(p_turb, calc_dB_some)
//       if (a){  // if True, calc_dB_some will point to the pure model
//          calc_dB_some(pos);
//       }
//       else{
//          calc_dB_SLAB(pos);
//          calc_dB_2D(pos);
//       }
//       --------------
//       De esta forma, no tenemos q chequear en c/iteracion!!
//       Tmb probar "const Doub const *pos".
// NOTA: 'pos' es el vector posicion (x,y,z) en [cm]
//
// TODO: Sino sacar los if( ... > 0.0 ) xq nunca los uso!
void MODEL_TURB::calc_dB(const Doub *pos){
    if(p_turb.percent_slab > 0.0) // solo calculo si vale la pena
        calc_dB_SLAB(pos);
    if(p_turb.percent_2d > 0.0)   // solo calculo si vale la pena
        calc_dB_2D(pos);

    for(int i=0; i<2; i++){       // solo las componentes "x, y" son turbulentas
        dB[i]   = dB_SLAB[i] + dB_2D[i];    // [G]
        //printf(" dB_2D[%d]:   %g\n", i, dB_2D[i]);
        //printf(" dB_SLAB[%d]: %g\n", i, dB_SLAB[i]); getchar();
    }
}


// TODO: probar "const Doub const *pos".
void MODEL_TURB::calc_B(const Doub *pos){
    calc_dB(pos);
    B[0] = dB[0];               // [G] Bx
    B[1] = dB[1];               // [G] By
    B[2] = dB[2] + p_turb.Bo;   // [G] Bz
}

#endif // DEFS_TURB_CC
/*------------------------------------------------------*/
