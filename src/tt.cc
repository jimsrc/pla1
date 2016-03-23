#ifndef TT_CC
#define TT_CC
#include "tt.h"

jim::jim(){
    a = 77;
    h1 = 0.077;
    aa = 0.9994;
}


int run(double x1, double rigidity){
    double aux=2.0;
    scl.build(rigidity);
    printf(" scl.vel: %g\n", scl.vel);

    return int(x1*aux);
}

#endif // TT_CC
