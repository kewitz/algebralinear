#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

#define EPSILON 1E-16

inline void assert(int test, char* msg) {
    if (!test) {
        printf(msg);
        exit(1);
    }
}

inline void SwapMLine(int n, int line1, int line2, double* A) {
    int j;
    double temp;
    for (j = 0; j < n; j++) {
        temp = A[line1*n + j];
        A[line1*n + j] = A[line2*n + j];
        A[line2*n + j] = temp;
    }

    return;
}

inline void SwapVLine(int index1, int index2, double* v) {
    double temp = v[index1];
    v[index1] = v[index2];
    v[index2] = temp;
}

/*
 * Método de Crout.
 */
void LUCroutDecompose(int n, double* A, double* L, double *U) {
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

void LUCroutInplaceDecompose(int n, double* A) {
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

double GaussPivot(int line, int n, double* A, double *b) {
    int i, j, k, iswap = 0;

    double temp, max = 0.0;
    for (i = line; i < n; i ++) {
        temp = abs(A[i*n + line]);
        if (temp > max) {
            max = temp;
            iswap = i;
        }
    }
    if (iswap > line) {
        SwapMLine(n, line, iswap, A);
        SwapVLine(line, iswap, b);
    }

    return A[line*n + line];
}

void GaussJordan(int n, double* A, double* x, double* b) {
    int i, j, k;
    double pivot, w;

    // Triangula.
    for (i = 0; i < n-1; i++) {
        pivot = GaussPivot(i, n, A, b);
        //assert(abs(pivot) > EPSILON, "Matriz singular.");

        for (j = i+1; j < n; j++) {
            w = A[j*n + i] / pivot;
            for (k = i+1; k < n; k++) {
                A[j*n + k] -= w * A[i*n + k];
            }
            b[j] -= w * b[i];
        }
    }

    // Retrosubstituição
    for (i = n-1; i >= 0; i--) {
        //assert(abs(pivot) > EPSILON, "Matriz singular.");
        x[i] = b[i];
        for (j = i+1; j < n; j++) {
            x[i] -= A[i*n + j] * x[j];
        }
        x[i] /= A[i*n + i];
    }
    return;
}
