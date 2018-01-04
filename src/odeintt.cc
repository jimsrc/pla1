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
    //WATCH_MEMORY 
    out.step_save = MatDoub(2,MAXSTP,0.0);
    #endif //MONIT_STEP

    #ifdef GUIDING_CENTER
    //WATCH_MEMORY 
    out.gc = new GuidingCenter(MAXSTP);
    #endif //GUIDING_CENTER
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
        // TODO: not exactly 'mu_old' here, review!
		out.out(-1,x,y,s,h,mu_old);				// aqui solo guarda x,y
	else{
		out.save(x,y);
		i++;
		cout << " i " << i << endl;}
	dtau = 0.0;

    #ifdef GUIDING_CENTER
    out.gc->n = 0; // total number of steps, recorded in 'out.gc'
    #endif //GUIDING_CENTER

	for (nstp=0;nstp<MAXSTP;nstp++) {
        // save current state
		save_history();					//--- scattering stuff

		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
        // advance one step, proposing the time-step 'h'
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
			out.out(nstp, x, y, s, s.hdid, mu_new);		// guarda solo si x>xout
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

    //--- ERROR: exit program.
    printf("\n ---> Too many steps in routine Odeint <--\n");
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
    
    #ifdef GUIDING_CENTER
    out.gc->calc_gc(&dydx[0], &y[0], x); // calculo cto de giro
    #endif //GUIDING_CENTER

	if((mu_old*mu_new)<0.0){
		out.nreb++;
		if(out.nreb>=out.nfilTau) 
            out.resizeTau();

		//--- guardo cosas de la "colisiones" con las irregularidades:
		out.Tau[out.nreb-1][0] = x; //[1] time @ collision
		out.Tau[out.nreb-1][1] = dtau; //[1] tiempo-de-colision instantaneo
		out.Tau[out.nreb-1][2] = NORM(y[0],y[2],0.0); //[1] posic "perpend" en q ocurre dicha "colision"
		out.Tau[out.nreb-1][3] = y[4]; //[1] posic "parall" en q ocurre dicha "colision"
        out.Tau[out.nreb-1][4] = acos(par.B[2]/Bmod)*180./M_PI; // [deg] angulo entre el vector 'B' y z. Siendo 0.0 para B-vector paralelo a versor positivo ^z+.

        //--- append this trail to 'ptrails'
        #ifdef WATCH_TRAIL
        // NOTE: before resetting 'dtau', let's see if it's close to
        // a resonance!
        //dtau_res = fmodf(dtau,2.*M_PI)/(2.*M_PI);
        // NOTE: 'dtau_res' will always be in the interval (0.0, 1.0), so 
        // in order to detect a 'bounce', we check whether 'dtau_res' is 
        // close to 0.0 or 1.0. And by "close" we mean it should be higher 
        // than 0.8, or lower that 0.1.
        // This means that the following condition is 'true' if dtau is 
        // true if 'dtau_res' is close to an integer multiple of 2*M_PI; i.e.
        // is true for all gryo-resonances!
        // The 'dtau>0.8*(2.*M_PI)' is to make sure that we are grabbing real
        // resonances. It says that the collision-time should be at least greater
        // than one gyrocycle (i.e. greater-ish than 2*M_PI).
        //if ( dtau>0.8*(2.*M_PI) && (dtau_res>0.8 || dtau_res<0.1) ){
        bool cond;
        Doub dtau_ = dtau/(2.*M_PI);    // collision-time normalized
        cond = false;
        for(int i=0; i<out.nbands; i++)
            cond |= (dtau_ > out.tau_bd[2*i]) && (dtau_ < out.tau_bd[2*i+1]);

        if (cond){
            out.append_trail(dtau); // append the trail and collision-time
        }
        #endif //WATCH_TRAIL

        // reset the collision-time
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
