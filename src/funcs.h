#ifndef FUNCS_H
#define FUNCS_H
//#include "control.h"
//#include "general.h"
#include "nr3.h"
#include "defs_turb.h"

//---------------------------------------------------
class PARAMS : public MODEL_TURB{
	public:
        #ifdef CYTHON
        PARAMS(){};// good to have for cython handling
        #endif //CYTHON
		PARAMS(string);
		//void calc_Bfield(VecDoub_I &);
		//PARAMS & operator=(const PARAMS &rhs);
	/*private:
		double pos[3];*/
};


#ifdef WATCH_TRAIL
class trail{
    public:
        trail();
        trail(int n, Doub tsize);
        void insert(const Doub *pos);
        ~trail();
        Doub **buffer;
        int n;
        Doub tsize;
        Doub dt;
};
#endif //WATCH_TRAIL


//---------------------------------------------------
struct rhs{  
    //functor for ode; copied from Numerical Recipes
    //      Doub eps;
    //      rhs(Doub epss) : eps(epss){}
    Doub bx, by, bz;
    void operator() (PARAMS par, const Doub x, VecDoub_I &y, VecDoub_O &dydx);
};


#ifdef GUIDING_CENTER
class GuidingCenter{
    public:
        GuidingCenter(Int len);
        void calc_gc(Doub* dydx, Doub* y, Doub x);
        MatDoub r_gc; // guiding center pos
        Doub* t; // time [1]
        Int n; // total number of steps
};
#endif //GUIDING_CENTER


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
		void build(const string, Int, Doub, Int, Int, int, int, char*);
		Output(void);
		//Output(string, const Int, char*); // mal implementado
		void init(const Int, const Doub, const Doub);
		void resize(void);
		void save_dense(Stepper &, const Doub, const Doub);
		void save(const Doub, VecDoub_I &);
		void out(const Int,const Doub,VecDoub_I &,Stepper &,const Doub);
		void save2file(void);
        void claim_own(void);
		bool file_exist(void);
		void resizeTau(void);

        #ifdef WATCH_TRAIL
        trail ptrail;
        Mat3DDoub ptrails;  // 
        Doub xtrail;        // stop-times to insert positions on trail
        void append_trail(); // append 'ptrail' to 'ptrails'
        Int ntrails;        // total number of appended trails
        #endif //WATCH_TRAIL

        #ifdef MONIT_SCATTERING
		//esto lo agrego para guardar cosas de la historia de
		//las trayectorias:
		int nfilTau, ncolTau;	// tamanio para 'Tau'
		int nreb;	 // nro de rebotes/scatterings en pitch
		MatDoub Tau; // tiempo de camino libre medio paralelo, y su posic x
		VecDoub mu;
        #endif //MONIT_SCATTERING

        #ifdef GUIDING_CENTER
        //MatDoub r_gc;
        GuidingCenter* gc;
        #endif //GUIDING_CENTER

		void set_Bmodel(PARAMS*);	// para apuntar al modelo q uso en main()
		void tic(void), toc(void);	// cronometro para c/pla
		Doub trun;			// tiempo de simulacion de c/pla
		Int nsteps;			// nro total de pasos de c/pla 
        // nombres de archivos de salida 
        char fname_out[200];
		char fname_trj[200];
        char fname_misc[200];
        char fname_owned[200];

        #ifdef MONIT_STEP
        MatDoub HistStep;
        MatDoub HistSeq;
        #if __cplusplus <= 199711L
        static const Doub MaxStep=1.0;
        #else
        static constexpr Doub MaxStep=1.0f;
        #endif //__cplusplus
        static const Int NStep=500;
        Doub dstep, dstep_part;
        //void monit_step(const Doub hdid);
        void monit_step(const Stepper s);
        void build_HistSeq(const Stepper s);
        MatDoub step_save;
        #endif //MONIT_STEP
        
	private:
		PARAMS *pm;
		Doub vmod, bmod;
		void save_pitch(void);
		Doub bx, by, bz, vx, vy, vz;
		VecDoub XSaveGen;	// tiempos de la salida
		Int n_tscales, cc, ndec, inid, nd, base, maxd;	
		string str_tscale; //tipo de escala temporal para la salida
		Doub decade, dt;
		void set_savetimes(Doub);
		ofstream ofile_trj;
        ofstream ofile_misc;
        ofstream ofile_own; 
    
    PRIVATE_OR_PUBLIC: // depends on CYTHON macro
		//----- histo del 'Tau'
		MatDoub HistTau;
		void build_HistTau(void);
		Doub dTau, maxTau, avrTau;
		Int nHistTau, nTau, dimHistTau;

        //--- histo del theta_coll
        Int nThColl; // has to be even!
        MatDoub HistThColl;
        void build_ThetaColl(void);
};


/*----- FUNCIONES NORMALES -----*/

double calc_gamma(double v);
void read_params(string fname, Doub &RIGIDITY, Doub &FRAC_GYROPERIOD, 
        Doub &NMAX_GYROPERIODS, Int &NPOINTS, Doub &ATOL, Doub &RTOL, 
        Int &n_Brealiz, string& str_timescale, Doub& tmaxHistTau, Int& nHist, Int& nThColl);
void init_orientation(int i, Doub **array_ori, VecDoub &y);
void LiberaMat(Doub **Mat, int i);
Doub **AllocMat(int nFilas, int nColumnas);
Doub **read_orientations(string fname, int &n);

#endif //FUNCS_H
//EOF
