#ifndef STEPPERBS_H
#define STEPPERBS_H
//#include "stepper.h"
#include "nr3.h"
#include "funcs.h"
//extern class PARAMS;

//----------------------------------------- struct StepperBase
struct StepperBase {
	Doub &x;
	Doub xold;
	VecDoub &y,&dydx;
	Doub atol,rtol;		// absolute, and relative tolerance
	bool dense;
	Doub hdid;
	Doub hnext;
	Doub EPS;
	Int n,neqn;
	VecDoub yout,yerr;
	StepperBase(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx, const Doub atoll,
		const Doub rtoll, bool dens);
};




//----------------------------------------- class StepperBS
template <class D>
struct StepperBS : StepperBase {
	typedef D Dtype;
	static const Int KMAXX=24;	// original:8
	static const Int IMAXX=KMAXX+1;
	Int k_targ;
	VecInt nseq;
	VecInt cost;
	MatDoub table;
	VecDoub dydxnew;
	Int mu;
	MatDoub coeff;
	VecDoub errfac;
	MatDoub ysave;
	MatDoub fsave;
	VecInt ipoint;
	VecDoub dens;
	StepperBS(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx, const Doub atol,
		const Doub rtol, bool dens, PARAMS parr);
	void step(const Doub htry,D &derivs);
	virtual void dy(VecDoub_I &y, const Doub htot, const Int k, VecDoub_O &yend,
		Int &ipt, D &derivs);
	void polyextr(const Int k, MatDoub_IO &table, VecDoub_IO &last);
	virtual void prepare_dense(const Doub h,VecDoub_I &dydxnew, VecDoub_I &ysav,
		VecDoub_I &scale, const Int k, Doub &error);
	virtual Doub dense_out(const Int i,const Doub x,const Doub h);
	virtual void dense_interp(const Int n, VecDoub_IO &y, const Int imit);
	PARAMS par;
};

#endif // STEPPERBS_H
//EOF
