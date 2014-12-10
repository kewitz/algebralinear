#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

extern void LUDec(int n, double* A, double* L, double* U) {
    LUCroutDecompose(n , A, L, U);
    return;
}
extern void LUInDec(int n, double* A) {
    LUCroutInplaceDecompose(n, A);
    return;
}
extern void SolveGJ(int n, double* A, double* x, double* b) {
    GaussJordan(n, A, x, b);
    return;
}
