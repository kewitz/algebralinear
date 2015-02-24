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

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#define DEBUG true

__global__ void JacobiIter(int n, int k, double* A, double* x, double* b) {
    int i = blockIdx.x * blockDim.x + threadIdx.x, j;
    if (i >= n) return;

    int off0 = ((k % 2) * n);
    int off1 = n - off0;

    double w = b[i];
    for (j = 0; j < i; j++) {
        w -= A[i*n + j] * x[j + off0];
    }
    for (j = i+1; j < n; j++) {
        w -= A[i*n + j] * x[j + off0];
    }
    w /= A[i*n + i];
    x[i + off1] = w;
}

extern "C" void CUJacobi(int n, int ks, double* A, double* x, double* b) {
    int k;

    double *dA, *dx, *db;
    if (DEBUG) printf("[+] CUDA Malloc...\n");
    cudaMalloc(&dA, sizeof(double)*n*n);
    cudaMalloc(&dx, sizeof(double)*n*2);
    cudaMalloc(&db, sizeof(double)*n);

    if (DEBUG) printf("[+] Copying to device memory...\n");
    cudaMemcpy(dA, A, sizeof(double)*n*n, cudaMemcpyHostToDevice);
    cudaMemcpy(dx, x, sizeof(double)*n, cudaMemcpyHostToDevice);
    cudaMemcpy(db, b, sizeof(double)*n, cudaMemcpyHostToDevice);

    const dim3 threads(64, 1);
    const dim3 blocks(1 + n/64, 1);
    if (DEBUG) printf("[+] Running Kernel...\n");
    for (k = 0; k < ks; k++) {
        JacobiIter<<<blocks, threads>>>(n, k, dA, dx, db);
        cudaDeviceSynchronize();
    }
    if (DEBUG) printf("[+] Copying Result and freeing memory...\n");
    cudaMemcpy(x, dx, sizeof(double)*n, cudaMemcpyDeviceToHost);

    cudaFree(dA);
    cudaFree(dx);
    cudaFree(db);
    if (DEBUG) printf("[+] Done.\n");
    return;
}
