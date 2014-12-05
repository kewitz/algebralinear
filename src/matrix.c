#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 1E-16

/*
 * Método de Crout.
 */
extern void LUCroutDecompose(int n, double* A, double* L, double *U) {
    int i, j, k;
    double soma;

    for (i = 0; i < n; i++) {
        U[i*n + i] = 1.0;
    }

    for (j = 0; j < n; j++) {
        for (i = j; i < n; i++) {
            soma = A[i*n + j];
            for (k = 0; k < j; k++) {
                soma -= L[i*n + k] * U[k*n + j];
            }
            L[i*n + j] = soma;
        }
        for (i = j; i < n; i++) {
            soma = 0.0;
            for (k = 0; k < j; k++) {
                soma += L[j*n + k] * U[k*n + i];
            }
            U[j*n + i] = (A[j*n + i] - soma) / L[j*n + j];
        }
    }
    return;
}

extern void LUCroutInplaceDecompose(int n, double* A) {
    int i, j, k;
    double soma, w;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            soma = A[i*n + j];
            for (k = 0; k < fmin(i, j); k ++) {
                soma -= A[i*n + k] * A[k*n + j];
            }
            if (j > i) {
                w = A[i*n + i];
                soma /= w;
            }
            A[i*n + j] = soma;
        }
    }
}

extern void GaussJordan(int n, double* A, double* b) {
    
}
