#ifndef FUNCS_H
#define FUNCS_H

//#include "control.h"
//#include "general.h"
#include "nr3.h"
#include "defs_turb.h"

//---------------------------------------------------
class PARAMS : public MODEL_TURB{
	public:
		PARAMS(string);
		void calc_Bfield(VecDoub_I &);
		//PARAMS & operator=(const PARAMS &rhs);
	private:
		double pos[3];
};



//---------------------------------------------------
struct rhs{  
    //functor for ode; copied from Numerical Recipes
    //      Doub eps;
    //      rhs(Doub epss) : eps(epss){}
    //void operator() (PARAMS par, const VecDoub x, VecDoub_I &y, VecDoub_O &dydx ){ 
    Doub bx, by, bz;
    void operator() (PARAMS par, const Doub x, VecDoub_I &y, VecDoub_O &dydx);
};



//---------------------------------------------------
template <class Stepper>
class Output {
	public:
		Int kmax;
		Int nvar;
		Int nsave;
		bool dense;
		Int count;
		Doub x1,x2,xout,dxout;
		VecDoub xsave;
		MatDoub ysave;
		//void build(string, Int, Doub, Int, char*); 
		void build(const string, Int, Doub, Int, int, int, char*);
		Output(void);
		//Output(string, const Int, char*); // mal implementado
		void init(const Int, const Doub, const Doub);
		void resize(void);
		void save_dense(Stepper &, const Doub, const Doub);
		void save(const Doub, VecDoub_I &);
		void out(const Int,const Doub,VecDoub_I &,Stepper &,const Doub);
		void save2file(void);
		ofstream ofile_trj;
        ofstream ofile_misc;
        ofstream ofile_own; 
        void claim_own(void);
		bool file_exist(void);
		void resizeTau(void);
		//esto lo agrego para guardar cosas de la historia de 
		//las trayectorias:
		int nfilTau, ncolTau;		// tamanio para 'Tau'
		int nreb;			// nro de rebotes/scatterings en pitch
		MatDoub Tau;			// tiempo de camino libre medio paralelo, y su posic x
		VecDoub mu;
		void set_Bmodel(PARAMS);	// para apuntar al modelo q uso en main()
		void tic(void), toc(void);	// cronometro para c/pla
		Doub trun;			// tiempo de simulacion de c/pla
		Int nsteps;			// nro total de pasos de c/pla 
        // nombres de archivos de salida 
        char fname_out[200];
		char fname_trj[200];
        char fname_misc[200];
        char fname_owned[200];
	private:
		PARAMS *pm;
		double pos[3], vmod, bmod;
		void save_pitch(void);
		double bx, by, bz, vx, vy, vz;
		VecDoub XSaveGen;	// tiempos de la salida
		int n_tscales, cc, ndec, inid, nd, base, maxd;	
		string str_tscale; //tipo de escala temporal para la salida
		double decade, dt;
		void set_savetimes(Doub);
		//----- histo del 'Tau'
		MatDoub HistTau;
		double dTau, maxTau, avrTau;
		int nHistTau, nTau, dimHistTau;
		void build_HistTau(void);
};


/*----- FUNCIONES NORMALES -----*/

double calc_gamma(double v);
void read_params(string fname, Doub &RIGIDITY, Doub &FRAC_GYROPERIOD, 
        Doub &NMAX_GYROPERIODS, int &NPOINTS, Doub &ATOL, Doub &RTOL, 
        int &n_Brealiz, string& str_timescale, Doub& tmaxHistTau, int& nHist);
void init_orientation(int i, Doub **array_ori, VecDoub &y);
void LiberaMat(Doub **Mat, int i);
Doub **AllocMat(int nFilas, int nColumnas);
Doub **read_orientations(string fname, int &n);

#endif //FUNCS_H
//EOF
