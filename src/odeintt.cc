#ifndef ODEINT_CC
#define ODEINT_CC
#include "nr3.h"
#include "funcs.h"
#include "odeintt.h"

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
	const Doub atol, const Doub rtol, const Doub h1, const Doub hminn,
	Output<Stepper> &outt,typename Stepper::Dtype &derivss, PARAMS parr, int wrankk) : 
	nvar(ystartt.size()),
	y(nvar),yold(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,rtol,dense, parr), 
	par(parr), wrank(wrankk) {
	EPS=numeric_limits<Doub>::epsilon();
	h=SIGN(h1,x2-x1);
	for (Int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn, x1, x2);

    #ifdef MONIT_STEP
    out.step_save = MatDoub(2,MAXSTP,0.0);
    #endif //MONIT_STEP

    #ifdef MONIT_SCATTERING
    out.gc = new GuidingCenter(MAXSTP);
    #endif
}


#ifdef KILL_HANDLER
template<class Stepper>
void Odeint<Stepper>::abort_mission(int signum){
    printf(" [r:%d] ---> ABORTANDO SIMULACION (signal: %d): %s\n", wrank, signum, out.fname_owned);
    char syscommand[4000];
    sprintf(syscommand, "rm %s", out.fname_owned);
    if (system(syscommand)==0)
        printf(" [r:%d] removed: %s\n", wrank, out.fname_owned);
    else
        printf(" [r:%d] Couldn't remove!: %s\n", wrank, out.fname_owned);
}
#endif //KILL_HANDLER


template<class Stepper>
void Odeint<Stepper>::integrate() {
    #ifdef MONIT_STEP
    //out.build_HistSeq(s);       // inicializa histog del "nseq"
    #endif

	int i=0;
	derivs(par, x, y, dydx);
	if (dense)
		out.out(-1,x,y,s,h);				// aqui solo guarda x,y
	else{
		out.save(x,y);
		i++;
		cout << " i " << i << endl;}
	dtau = 0.0;
    out.gc->n = 0; // total number of steps, recorded in 'out.gc'
	for (nstp=0;nstp<MAXSTP;nstp++) {
		save_history();					//--- scattering stuff
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h, derivs);

        #ifdef MONIT_SCATTERING
		check_scattering();				//--- scattering stuff
        #endif //MONIT_SCATTERING

		if (s.hdid == h) ++nok; else ++nbad;

        #ifdef MONIT_STEP
        out.step_save[0][nstp] = s.hdid;            // total step-size
        out.step_save[1][nstp] = s.hdid/s.nstep;    // partial (true) step-size
        #endif //MONIT_STEP

		if (dense){
			out.out(nstp, x, y, s, s.hdid);		// guarda solo si x>xout
		}
		else
			out.save(x,y);

		if ((x-x2)*(x2-x1) >= 0.0) {			// aqui termina la simulacion
			for (Int i=0;i<nvar;i++) ystart[i]=y[i];
			if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS){
				out.save(x,y);
			}
			out.nsteps = nstp;			// me gusta saber el nro total de pasos
			return;
		}
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
	throw("Too many steps in routine Odeint");
}


#ifdef MONIT_SCATTERING
template<class Stepper>
void Odeint<Stepper>::check_scattering(){
    Doub xyz[3] = {y[0],y[2],y[4]};
	par.calc_B(xyz);
	Bmod = NORM(par.B[0],par.B[1],par.B[2]);
	vmod = NORM(y[1],y[3],y[5]); 
	mu_new = (y[1]*par.B[0] + y[3]*par.B[1] + y[5]*par.B[2])/(vmod*Bmod);
	//-------------------------
	dtau += s.hdid;		// controlo cuanto pasa hasta el prox rebote
	//-------------------------
    out.gc->calc_gc(&dydx[0], &y[0], x); // calculo cto de giro
	if((mu_old*mu_new)<0.0){
		out.nreb++;
		if(out.nreb>=out.nfilTau) 
            out.resizeTau();

		// guardo cosas de la "colisiones" con las irregularidades:
		out.Tau[out.nreb-1][0] = x; //[1] time @ collision
		out.Tau[out.nreb-1][1] = dtau; //[1] tiempo-de-colision instantaneo
		out.Tau[out.nreb-1][2] = NORM(y[0],y[2],0.0); //[1] posic "perpend" en q ocurre dicha "colision"
		out.Tau[out.nreb-1][3] = y[4]; //[1] posic "parall" en q ocurre dicha "colision"
        out.Tau[out.nreb-1][4] = acos(par.B[2]/Bmod)*180./M_PI; // [deg] angulo entre plano x-y y z. Siendo 0.0 para B-vector paralelo a versor positivo ^z+.
		dtau = 0.0;
	}
}
#endif //MONIT_SCATTERING


template<class Stepper>
void Odeint<Stepper>::save_history(){
    Doub xyz[3] = {y[0],y[2],y[4]};
	par.calc_B(xyz);
	Bmod = NORM(par.B[0],par.B[1],par.B[2]);	// [1]
	vmod = NORM(y[1],y[3],y[5]); // [1]
	mu_old = (y[1]*par.B[0] + y[3]*par.B[1] + y[5]*par.B[2])/(vmod*Bmod);
}


// declare/define (?) class w/ specific template (BS: Bulirsch-Stoer)
#include "stepperbs.h"
template class Odeint<StepperBS<rhs> >;


#endif // ODEINT_CC
//EOF
