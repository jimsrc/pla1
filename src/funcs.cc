#ifndef FUNCS_CC
#define FUNCS_CC

#include "control.h"
#include "funcs.h"
#include "general.h"

using namespace std;


// declaration of definition in general.cc
//extern ESCALAS scl;


/*----- FUNCIONES NORMALES -----*/

// recibe una velocidad adimensionalizada
// TODO: convertir esto en inline o macro!
double calc_gamma(double v){
	double beta, gamma;
	beta = v*scl.vel / clight;
    #ifdef BETA_CHECK
    if (beta>=1.0)
        printf(" beta>=1.0!!, beta: %g, v: %g\n", beta, v);
    #endif
	gamma = pow(1. - beta*beta, -.5);
	return gamma;
}


/* lee parametros input en main() */
void read_params(string fname, Doub &RIGIDITY, Doub &FRAC_GYROPERIOD, 
        Doub &NMAX_GYROPERIODS, Int &NPOINTS, Doub &ATOL, Doub &RTOL, 
        Int &n_Brealiz, string& str_timescale, Doub& tmaxHistTau, 
        Int& nHist, Int& nThColl){
	string dummy;
	ifstream filein(fname.c_str());
	if (!filein.good()) {
		cout << " problema al abrir " << fname << endl;
		exit(1);
	}
	filein >> RIGIDITY		>> dummy;  // [V] rigidez de las plas
	filein >> FRAC_GYROPERIOD	>> dummy;  // [1] fraction of gyroper
	filein >> NMAX_GYROPERIODS	>> dummy;  // [1] nro of gyroperiods
	filein >> NPOINTS		>> dummy;  // [1] nro output pts
	filein >> ATOL			>> dummy;  // [units of *y] abs tolerance
	filein >> RTOL			>> dummy;  // [1] rel tolerance
	filein >> n_Brealiz		>> dummy;  // [1] nmbr of B-field realizations
	filein >> str_timescale >> dummy;  // [string] nombre de la escala temporal para la salida
	filein >> tmaxHistTau	>> dummy;  // [1] max collision time (in gyro-period units) for histogram of collision times 
	filein >> nHist			>> dummy;  // [1] nmbr of bins for histogram of collision-times
	filein >> nThColl		>> dummy;  // [1] nmbr of bins for histogram of the angle between the x-y plane and z axis.
}



void init_orientation(int i, Doub **array_ori, VecDoub &y){
	double th, ph, mu;	// theta, phi y pitch-cosine
	th	= array_ori[i][0];
	ph	= array_ori[i][1];
	mu	= cos(th);	// pitch

	y[1]	= sqrt(1.-mu*mu)*cos(ph);	// [1] vx
	y[3]	= sqrt(1.-mu*mu)*sin(ph);	// [1] vy
	y[5]	= mu;				// [1] vz
	// en el origen siempre
	y[0]	= 0.0;				// [1] x
	y[2]	= 0.0;				// [1] y
	y[4]	= 0.0;				// [1] z
}



void LiberaMat(Doub **Mat, int i){ 
    for(int k=0; k<=i; k++)
            free(Mat[i]);
    free(Mat);
}



Doub **AllocMat(int nFilas, int nColumnas){
    Doub **Mat;

    if((Mat = (Doub**) calloc(nFilas, sizeof(Doub *))) == NULL)
        return NULL;

    for(int i=0; i<nFilas; i++)
        if((Mat[i] = (Doub *) calloc(nColumnas, sizeof(Doub)))==NULL) {
            LiberaMat(Mat,i);
            return NULL;
        }   
    return Mat;
}



Doub **read_orientations(string fname, int &n){
	double dummy;
	ifstream filein(fname.c_str());
	if (!filein.good()) {
		cout << " problema al abrir " << fname << endl;
		exit(1);
	}
	n = 0;
	for(;;n++){
		filein >> dummy	>> dummy;
		if(filein.eof()) break;
	}

	ifstream file(fname.c_str());

	double **array;
	array	= AllocMat(n, 2);
	for(int i=0; i<n; i++){
		file >> array[i][0] >> array[i][1];
	}
	return array;
}



//----------------------- class Ouput
//void Output<Stepper>::build(string str_tscalee, Int nsavee, Doub tmaxHistTau, Int nHist, char* fname_out){ 
template <class Stepper>
void Output<Stepper>::build(const string str_tscalee, Int nsavee, Doub tmaxHistTau, Int nHist, Int nThColl_, int i, int j, char *dir_out){
	kmax	= 500;
	nsave	= nsavee;
	count	= 0;
	xsave.resize(kmax);
	dense 	= nsave > 0 ? true : false;
	//-------------------- archivos de salida 
    sprintf(fname_out, "%s/B%02d_pla%03d", dir_out, j, i);
	sprintf(fname_trj,  "%s_traj.dat",  fname_out);
	sprintf(fname_misc, "%s_misc.dat", fname_out);
	sprintf(fname_owned, "%s.owned", fname_out);

	str_tscale = str_tscalee;	// tipo de escala temporal para la salida

    #ifdef MONIT_SCATTERING
	//-------- cosas de scatterings:
	nfilTau = 500;
	nreb	= 0;  // nro inic de rebotes
	ncolTau	= 4;  // 4 columnas: 1 para el tau de scattering, 2 para las posic parall/perp, y 1 para el angulo entre el plano x-y y z.
	Tau 	= MatDoub(nfilTau, ncolTau, 0.0); // (*) para grabbar tiempos de scattering, y la posic x
	// (*): inicializo en ceros

	//-------- histograma del 'Tau'
	nHistTau 	= nHist;				// nro de bines
	dTau 		= tmaxHistTau / nHistTau;		// ancho del bin
	dimHistTau	= 2;
	HistTau		= MatDoub(nHistTau, dimHistTau, 0.0);	// histog 1-D
	nsteps		= 0;

    //-------- histograma del 'ThetaColl'
    if(fmod(nThColl_, 2.0)!=0.0) {
        printf("\n ---> ERROR: 'nThColl' has to be even!!\n");
        exit(1);
    }
    nThColl     = nThColl_;
    HistThColl  = MatDoub(nThColl, 2, 0.0); // histog 1-D
    #endif //MONIT_SCATTERING

    #ifdef MONIT_STEP
    HistStep    = MatDoub(NStep, 4);
    dstep       = MaxStep/(1.0*NStep);
    dstep_part  = dstep/8.;
    //MinStep     = 1e3;
    for(int i=0; i<NStep; i++){
        HistStep[i][0] = (i+.5)*dstep;
        HistStep[i][1] = 0.0;               // counts
        HistStep[i][2] = (i+.5)*(dstep_part);
        HistStep[i][3] = 0.0;               // counts
    }
    printf(" NStep: %d\n", NStep);
    printf(" MaxStep: %g\n", MaxStep);
    //printf(" HistStep: %g\n", HistStep[0][1]);
    #endif //MONIT_STEP
}


template <class Stepper>
Output<Stepper>::Output() : kmax(-1),dense(false),count(0) {}

/*
// TODO: arreglar esta implementacion si la vas a usar
template <class Stepper>
Output<Stepper>::Output(string str_tscalee, const Int nsavee, char* fname){
	build(str_tscalee, nsavee, fname);
	//kmax(500),nsave(nsavee),count(0),xsave(kmax) {
	//dense = nsave > 0 ? true : false;
}
*/

template <class Stepper>
void Output<Stepper>::set_savetimes(Doub xhi){
	if(str_tscale=="linear"){
		dxout=(x2-x1)/nsave;
		mu  	 = VecDoub(nsave, 0.0);		// (*) pitch en los tiempos "xsave"
		XSaveGen = VecDoub(nsave, 0.0);
		for(int i=0; i<nsave; i++){
			XSaveGen[i] = (i+1)*dxout;
			//printf(" XMix(%d) [wc-1]: %g\n", i, XSaveGen[i]);
		}
		cc = 0;					// indice para 'XSaveGen'
	}
	else if(str_tscale=="mixed"){
		inid = 1;
		maxd = int(M_LOG10E*log(xhi));
		ndec = maxd - inid + 1;
		mu  	 = VecDoub(nsave*ndec, 0.0);		// (*) pitch en los tiempos "xsave"
		XSaveGen = VecDoub(nsave*ndec, 0.0);

		nd = inid;
		for(nd=inid; nd<maxd; nd++){
			dt = (pow(10, (1.0*(nd+1))) - pow(10, (1.0*nd)))/nsave;
			for(int i=0; i<nsave; i++){
				cc = i+(nd-inid)*nsave;
				XSaveGen[cc] = pow(10, 1.0*(nd)) + (i+1)*dt;
				//printf(" XMix(%d) [wc-1]: %g\n", cc, XSaveGen[cc]);
			}
		}

		dt = (xhi - pow(10, 1.0*(maxd) ))/nsave;
		for(int i=0; i<nsave; i++){
			cc = i+(maxd-inid)*nsave;
			XSaveGen[cc] = pow(10, 1.0*(maxd) ) + (i+1)*dt;
			//printf(" XMix(%d) [wc-1]: %g\n", cc, XSaveGen[cc]);
		}
		// reseteo indice:
		cc = 0;
	}
	else
		throw(" USAR 'linear' O 'mixed' SOLAMENTE! (Jimmy)");
}

template <class Stepper>
void Output<Stepper>::init(const Int neqn, const Doub xlo, const Doub xhi) {
	nvar=neqn;
	if (kmax == -1) return;
	ysave.resize(nvar,kmax);
	if (dense) {
		x1=xlo;
		x2=xhi;
		xout=x1;
		if(xout>0.0) throw(" NO ESTA IMPLEMENTADO PARA EMPEZAR EN TIEMPO t>0 (Jimmy).");
		//---- seteo los tiempos a guardar
		set_savetimes(xhi);
	}
}

template <class Stepper>
void Output<Stepper>::resize(){	// redimensiona el vector 'xsave' hacia el doble de su longitud, y 
				// redimensiona el array 'ysave' hacia las dimensiones (nvar,kmax).
				// Como no preserva los valores, los guarda temporalmente antes de 
				// redimensionar.Despues de redimensionar,recupera la data temporal.
	Int kold=kmax;
	kmax *= 2;
	VecDoub tempvec(xsave);	// backup de 'xsave'
	xsave.resize(kmax);
	for (Int k=0; k<kold; k++)
		xsave[k]=tempvec[k];
	MatDoub tempmat(ysave);	// backup de 'ysave'
	ysave.resize(nvar,kmax);
	for (Int i=0; i<nvar; i++)
		for (Int k=0; k<kold; k++)
			ysave[i][k]=tempmat[i][k];
}

template <class Stepper>
void Output<Stepper>::resizeTau(){	// redimensiona el vector 'xsave' hacia el doble de su longitud.
				// Como no preserva los valores, los guarda temporalmente antes de 
				// redimensionar.Despues de redimensionar,recupera la data temporal.
	Int nold=nfilTau;
	nfilTau *= 2;
	MatDoub tempmat(Tau);	// backup de 'Tau'
	Tau.resize(nfilTau, ncolTau);
	for(Int i=0; i<nold; i++)
		for(Int j=0; j<ncolTau; j++)
			Tau[i][j] = tempmat[i][j];
}


#ifdef MONIT_STEP
template <class Stepper>
void Output<Stepper>::build_HistSeq(const Stepper s){
    HistSeq = MatDoub(s.IMAXX, 2);
    for(int i=0; i<s.IMAXX; i++){
        HistSeq[i][0] = s.nseq[i];
        HistSeq[i][1] = 0.0;
    }
}


template <class Stepper>
//void Output<Stepper>::monit_step(const Doub hdid){
void Output<Stepper>::monit_step(const Stepper s){
    /*if(hdid!=0.05)
        MinStep = MIN(MinStep, hdid);*/
    int ns;
    if(s.hdid<=MaxStep){
        // h total
        ns = int(s.hdid/dstep);
        HistStep[ns][1]++;
    }
    if((s.hdid/s.nstep)<=(NStep*dstep_part)){
        // h partial
        ns = int((s.hdid/s.nstep)/(dstep_part));
        HistStep[ns][3]++;
    }
    //
    for(int i=0; i<HistSeq.nrows(); i++){
        if(s.nstep==HistSeq[i][0])
            HistSeq[i][1]++;
    }
}
#endif //MONIT_STEP


template <class Stepper>
void Output<Stepper>::save_pitch(){
	for(int i=0;i<3;i++)
		pos[i] = (ysave[(2*i)][count]) *scl.rl;	// [cm]
    //printf(" >>> save_pitch...\n");
	pm->calc_B(pos);

	bx=pm->B[0];		by=pm->B[1];		bz=pm->B[2];		// [G]
	vx=ysave[1][count];	vy=ysave[3][count];	vz=ysave[5][count];	// [1]

	bmod = pow(bx*bx + by*by + bz*bz, .5);
	vmod = pow(vx*vx + vy*vy + vz*vz, .5);
	mu[count] = vx*bx + vy*by + vz*bz;
	mu[count] /= vmod*bmod;
}

template <class Stepper>
void Output<Stepper>::save_dense(Stepper &s, const Doub xout, const Doub h){
	if (count == kmax) resize();
	for (Int i=0;i<nvar;i++){
		ysave[i][count]=s.dense_out(i,xout,h);
	}
	save_pitch();			// calcula mu
	xsave[count++]=xout;		// <=> xsave[count]=xout; count++;
}

template <class Stepper>
void Output<Stepper>::save(const Doub x, VecDoub_I &y) {
	if (kmax <= 0) return;
	if (count == kmax) resize();
	for (Int i=0;i<nvar;i++)
		ysave[i][count]=y[i];
	save_pitch();			// calcula mu
	xsave[count++] = x;//x;		// <=> xsave[count]=x; count++;
	//cc++;
}

template <class Stepper>
void Output<Stepper>::out(const Int nstp,const Doub x,VecDoub_I &y,Stepper &s,const Doub h) {
	//if(count>=200) {printf(" COUNT=%d AQUI!\n", count); getchar();}
	if (!dense)
		throw("dense output not set in Output!");
	if (nstp == -1) {
		save(x, y); 
		xout = XSaveGen[cc];
		cc++;	//+= dxout;
	} else {
		while ((x-xout)*(x2-x1) > 0.0) {
			save_dense(s, xout, h);		// interpola a 'xout'
			xout = XSaveGen[cc]; //+= dxout;			// avanza 'xout' en 'dxout'
			cc++;
		}
	}
}


//esto lo agrego para guardar cosas de la historia de 
//las trayectorias:
template <class Stepper>
void Output<Stepper>::set_Bmodel(PARAMS *pmm){
	pm = pmm;
}


// escribimos archivo dummy "owned" para flagear de q YO 
// estoy trabajando con esta pla
template <class Stepper>
void Output<Stepper>::claim_own(){
    ofile_own.open(fname_owned);
    ofile_own << "dummy" << endl;
    ofile_own.close();
}


// chekea si existe (aunq su tamanio sea 0 bytes) o no un archivo.
template <class Stepper>
bool Output<Stepper>::file_exist(){
   	//if(ifstream(fname_trj)){
   	if(ifstream(fname_owned)){
		printf("\n YA EXISTE: %s\n", fname_owned);
		return true;
	}
	printf("\n AUN NO EXISTE: %s\n", fname_owned);
	return false;
}


template <class Stepper>
void Output<Stepper>::build_HistTau(){
	maxTau = dTau*nHistTau;

	for(int i=0; i<nHistTau; i++){
		HistTau[i][0] = (i+.5)*dTau;	// bines centrados
	}

	for(int i=0; i<nreb; i++){
		if(Tau[i][0] <= maxTau){
			nTau = int(Tau[i][0] / dTau);
			HistTau[nTau][1] += 1;
		}
	}
	//--- calculo la media
	avrTau = 0.0;
	for(int i=0; i<nreb; i++){
		avrTau += Tau[i][0];
	}
	avrTau /= nreb;
}


template <class Stepper>
void Output<Stepper>::build_ThetaColl(){
    /* Warning: 'nThColl' has to be even! */
    int nth;
    Doub dth = 180.0/nThColl; // resolucion del histo

    // dominio en th=(-90,90) [deg]
    for(int i=0; i<nThColl; i++){
        HistThColl[i][0] = -90.0 + (i+.5)*dth; // [deg]
    }

    // Tau[:][3] ---> theta de colision
    for(int i=0; i<nreb; i++){
        nth =  int(Tau[i][3]/dth);
        nth += nThColl/2; // correction to avoid negative indexes
        HistThColl[nth][1]++;
    }
}


template <class Stepper>
void Output<Stepper>::save2file(){
	double t, x, y, z, v;
	double vx, vy, vz;
    double err, gamma;
	//-------------------- guardo la trayectoria
	//printf(" COUNT @ save2file: %d\n", count); getchar();
	ofile_trj.open(fname_trj);
	for(int i=0; i<count; i++){
		t 	= xsave[i] / scl.wc;			// [seg]
		x 	= ysave[0][i] * scl.rl/AU_in_cm; 	// [AU]
		y 	= ysave[2][i] * scl.rl/AU_in_cm; 	// [AU]
		z 	= ysave[4][i] * scl.rl/AU_in_cm; 	// [AU]
		vx 	= ysave[1][i];  			// [1] 
		vy 	= ysave[3][i];  			// [1]
		vz 	= ysave[5][i];  			// [1]
		v 	= pow(vx*vx + vy*vy + vz*vz, .5); 	// [1] TODO: CAMBIARLO A 'sqrt'
		gamma	= calc_gamma(v);
		err	= gamma/scl.gamma - 1.;			// error relativ del gamma relativista

		ofile_trj << setiosflags(ios :: showpoint | ios :: uppercase);
		ofile_trj << setw(5) << setprecision(8) << t << " ";
		ofile_trj << setw(5) << setprecision(8) << x << " ";    
		ofile_trj << setw(5) << setprecision(8) << y << " ";   
		ofile_trj << setw(5) << setprecision(8) << z << " ";
		ofile_trj << setw(5) << setprecision(8) << mu[i] << " ";
		ofile_trj << setw(5) << setprecision(8) << err << endl;
	}
	// cerramos archivo de trayectoria
	ofile_trj.close();

	/**** guardamos otras cosas sobre la historia de la trayectoria ***/
	//--- nro de rebotes, y colission-time promedio
	ofile_misc.open(fname_misc);
    ofile_misc << "# ***** MISC INFO *****" << endl;

    #ifdef MONIT_SCATTERING
    ofile_misc << "# Histogram on measured collision-times 'Tau'" << endl;
	ofile_misc << "# nro_rebotes: " << nreb << endl;
	//--- histograma de 'Taus'
	build_HistTau();
	ofile_misc << "# average_Tau: " <<setw(10)<<setprecision(8)<< avrTau << endl;	// [1]
	ofile_misc << "# trun/min: " <<setw(10)<<setprecision(8)<< (trun/60.) << endl;	// [sec]
	ofile_misc << "# steps: " <<setw(10)<<setprecision(10)<< nsteps << endl;
	ofile_misc << "#####" << endl;	// cadena para separar tipos de dato q grabo
	for(int i=0; i<nHistTau; i++){
		ofile_misc << HistTau[i][0] << " ";			// [1] bin centrado
		ofile_misc << setw(10) << HistTau[i][1] << endl;	// [1] nro de cuentas en este bin
	}

    //--- histograma del theta-en-colision
    build_ThetaColl();
	ofile_misc << "#####" << endl;	// cadena para separar tipos de dato q grabo
    ofile_misc << "# Histogram on angle between x-y plane and z axis (Theta_Coll)"<< endl;
	ofile_misc << "#####" << endl;	// cadena para separar tipos de dato q grabo
    for(Int i=0; i<nThColl; i++){
        ofile_misc << HistThColl[i][0] << " "; // [deg] bin centrado
        ofile_misc << setw(10) << HistThColl[i][1] << endl; // [1] nro de cuentas
    }
    #else
    ofile_misc << "# +++++ NO SCATTERING INFORMATION +++++" << endl;
    #endif //MONIT_SCATTERING

	// cerramos archivo de misc
	ofile_misc.close();
}

template <class Stepper>
void Output<Stepper>::tic(){
	trun = time(NULL);
}

template <class Stepper>
void Output<Stepper>::toc(){
	trun = time(NULL) - trun; // [sec] tiempo de corrida para 1 pla
}



//------------------------------------------- class PARAMS
PARAMS::PARAMS(string fname_turb):
	MODEL_TURB(fname_turb) {
}


void PARAMS::calc_Bfield(VecDoub_I &y){
    //printf(" >>> calc_Bfield...\n");
	pos[0] = y[0] *scl.rl;		// [cm] x
	pos[1] = y[2] *scl.rl;		// [cm] y
	pos[2] = y[4] *scl.rl;		// [cm] z
	calc_B(pos);
}



//-------------------------------------------
void rhs::operator() (PARAMS par, const Doub x, VecDoub_I &y, VecDoub_O &dydx ){
    //double bx, by, bz; 
    par.calc_Bfield(y);

    bx = par.B[0] / scl.Bo;
    by = par.B[1] / scl.Bo;
    bz = par.B[2] / scl.Bo;
    // rewrite x^2y"(x)+xy'(x)+x^2y=0 as coupled FOODEa
    dydx[0] = y[1];
    dydx[1] = y[3] * bz - y[5] * by; 
    dydx[2] = y[3];
    dydx[3] =-y[1] * bz + y[5] * bx; 
    dydx[4] = y[5];
    dydx[5] =-y[3] * bx + y[1] * by; 
}


// declare/define (?) class w/ specific template
#include "stepperbs.h"
template class Output<StepperBS<rhs> >; // rhs: system of equations I use!


#endif //FUNCS_CC
//EOF
