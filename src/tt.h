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

    private:
        double h1, h2;
};

extern ESCALAS scl;
int run(double x1, double rigidity);

#endif // TT_H 
