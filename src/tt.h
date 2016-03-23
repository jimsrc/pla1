#ifndef TT_H
#define TT_H
#include "control.h"
#include "general.h"
#include "defs_turb.h"
#include "funcs.h"

class jim{
    public:
        jim();
        int a, b;
        double aa;
        PARAMS_SEM s;

    private:
        double h1, h2;
};

extern ESCALAS scl;
int run(double x1, double rigidity);
int run2(double x1, double rigidity, PARAMS_TURB* pt);
int run3(double x1, double rigidity, PARAMS_TURB pt);

#endif // TT_H 
