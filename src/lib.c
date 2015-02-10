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

extern void LUDec(int n, double* A, double* L, double* U) {
    LUCroutDecompose(n , A, L, U);
    return;
}
extern void LUInDec(int n, double* A) {
    LUCroutInplaceDecompose(n, A);
    return;
}
extern void LUSolve(int n, double* LU, double* x, double* b) {
    LUSubstitute(n, LU, x, b);
    return;
}

extern void SolveGJ(int n, double* A, double* x, double* b) {
    GaussJordan(n, A, x, b);
    return;
}

extern void SolveJacobi(int n, int ks, double* A, double* x, double* b) {
    Jacobi(n, ks, A, x, b);
    return;
}

extern void SolveGaussSeidel(int n, int ks, double* A, double* x, double* b) {
    GaussSeidel(n, ks, A, x, b);
    return;
}
