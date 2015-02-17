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
 *
 */

#include <cuda.h>
#include <cuda_runtime.h>
#include "device.h"

//  Utils
int getCUDAdevices() {
    int deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
    return deviceCount;
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