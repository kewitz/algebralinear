/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Leonardo Kewitz
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

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

void LUSubstitute(int n, double* LU, double* x, double* b) {
    int i, j;
    double w, s, pivot, _eps = 1E-16;

    for (i = 0; i < n; i++) {
        w = x[i];
        for (j = 0; j < i; j++) {
            w -= LU[i*n + j] * x[j];
        }
        pivot = LU[i*n + i];
        assert(fabs(w) > _eps, "Matriz singular.");
        w /= pivot;
        x[i] = w;
    }

    for (i = n - 1; i >= 0; i--) {
        s = x[i];
        for (j = i+1; j < n; j++) {
            s -=  LU[i*n + j] * x[j];
        }
        x[i] = s;
    }
    return;
}

double GaussPivot(int col, int n, double* A, double *b) {
    int i, j, k, iswap = 0;

    double temp, max = 0.0;
    for (i = col; i < n; i ++) {
        temp = abs(A[i*n + col]);
        if (temp > max) {
            max = temp;
            iswap = i;
        }
    }
    if (iswap > col) {
        SwapMLine(n, col, iswap, A);
        SwapVLine(col, iswap, b);
    }

    return A[iswap*n + iswap];
}

void GaussJordan(int n, double* A, double* x, double* b) {
    int i, j, k;
    double pivot, w, _eps = 1E-16;

    // Triangula.
    for (j = 0; j < n-1; j++) {
        pivot = GaussPivot(j, n, A, b);
        assert(fabs(pivot) > _eps, "Matriz singular.");
        for (i = j+1; i < n; i++) {
            w = A[i*n + j] / pivot;
            for (k = j+1; k < n; k++) {
                A[i*n + k] -= w * A[j*n + k];
            }
            b[i] -= w * b[j];
        }
    }

    // Retrosubstituição
    for (i = n-1; i >= 0; i--) {
        assert(fabs(pivot) > _eps, "Matriz singular.");
        x[i] = b[i];
        for (j = i+1; j < n; j++) {
            x[i] -= A[i*n + j] * x[j];
        }
        x[i] /= A[i*n + i];
    }
    return;
}

void Jacobi(int n, int ks, double* A, double* x, double* b) {
    int i, j, k;
    double w, pivot, _alpha = 1E-8, _eps = 1E-16;
    double *xnext;

    xnext = malloc(n*sizeof(double));

    for (k = 0; k < ks; k++) {
        for (i = 0; i < n; i++) {
            pivot = A[i*n + i];
            w = b[i];
            assert(fabs(pivot) > _eps, "Matriz singular.");

            for (j=0; j < i; j++) {
                w -= A[i*n + j] * x[j];
            }
            for (j=i+1; j < n; j++) {
                w -= A[i*n + j] * x[j];
            }

            xnext[i] = w/pivot;
        }

        for (i = 0; i < n; i++) {
            x[i] = xnext[i];
        }
    }
}
