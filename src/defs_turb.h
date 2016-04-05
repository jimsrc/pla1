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
		Int n_modos;

	public:
		FASES(void) {};
		//~FASES(void);
		void build(int, PARAMS_SEM *);
		Doub *phi_s, *a_s, *b_s;
		Doub *phi_2d, *b_2d;
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
		const static Doub gS=(5./3.);		// potencia espectral slab
		const static Doub g2D=(8./3.);		// potencia espectral 2D
	public:
		PARAMS_TURB(string);	// constructor
		PARAMS_TURB(void);			// constructor (otro)
		//~PARAMS_TURB(void);			// descrtructor

		string FNAME_INPUT;

		Int n_modos;
		Doub lambda_min;
        Doub lambda_max;
		Doub Lc_slab, Lc_2d;	// longitudes de correlacion
		Doub sigma_Bo_ratio;
		Doub percent_slab;
		Doub percent_2d;

		Doub Bo;
		Doub sigma_S;
		Doub sigma_2D;
		Doub *dk;
		Doub *k;
		Doub *Bk_SLAB;
		Doub *Bk_2D;

		PARAMS_SEM sem;
		FASES fases;

        void build_spectra();
		void build(string);			// puedo usarlo si es q use el contructor con "void"
};

/*------------------ parametros turbulencia -----------------------*/
class MODEL_TURB{
	/*private:
		PARAMS_TURB p_turb;*/
	public:
		string FNAME_INPUT;
		MODEL_TURB(string fname_input); //{build(fname_input);};  // constructor
		MODEL_TURB(void) {};			// constructor "trivial"
		//~MODEL_TURB(void);			// destructor

		void build(string);
		void calc_dB_SLAB(const Doub *);
		void calc_dB_2D(const Doub *);
		void calc_dB(const Doub *);
		void calc_B(const Doub *);

		PARAMS_TURB params_turb(void);

		Doub *B;		// [G]
		Doub *dB;		// [G]
		Doub *dB_SLAB;	// [G]
		Doub *dB_2D;		// [G]
		PARAMS_TURB p_turb;
		void next_B_realization(void);
        void fix_B_realization(const int); // fija la realizacion en funcion del argumento
};

#endif //DEFS_TURB_H
