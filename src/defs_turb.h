#ifndef DEFS_TURB_H
#define DEFS_TURB_H
#include "nr3.h"
#include "general.h"
/*-------------------------------------------------------*/

class PARAMS_SEM{
	public:
		PARAMS_SEM(void) {};
		~PARAMS_SEM(void) {};
		long slab[3];
		long two[2];
};

/*-------------------------- fases random -------------------------------*/
class FASES{
	private:
		Int Nm_slab;
        Int Nm_2d;

	public:
		FASES(void) {};
		//~FASES(void);
		void build(Int, Int, PARAMS_SEM *);
		Doub *phi_s, *a_s, *b_s; // fases for slab
		Doub *phi_2d, *b_2d;     // fases for 2d
		void construir_fases_random(PARAMS_SEM);
};

/*-------------------------- parametros turbulencia --------------------------------*/
class PARAMS_TURB{
	private:
		void read_params(string);
		void build_sigmas(void);
		void build_k_and_dk(void);
		void build_Bk_SLAB(void);
		void build_Bk_2D(void);
		void report(void);
	public:
		PARAMS_TURB(string);	// constructor
		PARAMS_TURB(void);			// constructor (otro)
		//~PARAMS_TURB(void);			// descrtructor

		string FNAME_INPUT;

		//Int n_modos;
		Int Nm_slab;
		Int Nm_2d;
        Doub lmin_s, lmax_s;    // escalas de turbulencia Slab
        Doub lmin_2d, lmax_2d;  // escalas de turbulencia 2D
		Doub Lc_slab, Lc_2d;	// longitudes de correlacion slab && 2d
		Doub sigma_Bo_ratio;
		Doub percent_slab;
		Doub percent_2d;

		//Doub Bo;
		Doub sigma_S;
		Doub sigma_2D;
		Doub *dk_s, *dk_2d;
		Doub *k_s, *k_2d;
		Doub *Bk_SLAB;
		Doub *Bk_2D;

		PARAMS_SEM sem;
		FASES fases;

        void build_spectra();
		void build(string);			// puedo usarlo si es q use el contructor con "void"
};

/*------------------ parametros turbulencia -----------------------*/
class MODEL_TURB{
	private:
		//PARAMS_TURB p_turb;
		void calc_dB_SLAB(const Doub *);
		void calc_dB_2D(const Doub *);
		void calc_dB(const Doub *);

	public:
		string FNAME_INPUT;
		MODEL_TURB(string fname_input); //{build(fname_input);}; //constructor
		MODEL_TURB(void) {};			// constructor "trivial"
		//~MODEL_TURB(void);			// destructor
		void build(string);
		void calc_B(const Doub *);
        void fix_B_realization(const int); // fija la realizacion en funcion del argumento

    PRIVATE_OR_PUBLIC:
		PARAMS_TURB p_turb; //public only for cython
		Doub *dB_SLAB;	// [G]
		Doub *dB_2D;		// [G]
		Doub *B;		// [G]
		Doub *dB;		// [G]
};

#endif //DEFS_TURB_H
