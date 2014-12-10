#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

void printMatrix(int, double*, char*);
void printVector(int, double*, char*);

int main() {
    int i, j;
    double w;
    double * A, * l, * u, * b, * x;

    printf("Alocando memória...\n");
    A = (double*) malloc(sizeof(double)*9);
    l = (double*) malloc(sizeof(double)*9);
    u = (double*) malloc(sizeof(double)*9);
    b = (double*) malloc(sizeof(double)*3);

    printf("Preenchendo matriz A...\n");
    for (i = 0; i < 9; i++) {
        w = (double) (rand()%100);
        A[i] = w;
    }
    printf("Preenchendo vetor b...\n");
    for (i = 0; i < 3; i++) {
        w = (double) (rand()%100);
        b[i] = w;
    }

    printf("Dados utilizados no teste:\n");
    printMatrix(3, A, "A");
    printVector(3, b, "b");

    printf("Executando teste de decomposição A=LU...\n");
    LUCroutDecompose(3, A, l, u);
    printMatrix(3, l, "L");
    printMatrix(3, u, "U");

    printf("Executando teste de eliminação de Gauss-Jordan...\n");
    x = (double*) malloc(sizeof(double)*3);
    GaussJordan(3, A, x, b);
    printMatrix(3, A, "A");
    printVector(3, x, "x");

    return 0;
}

void printVector(int size, double* data, char * name) {
    int i;
    printf("%s:", name);
    for (i = 0; i < size; i++) {
        printf("\t%.3f\n", data[i]);
    }
    printf("\n");
}
void printMatrix(int size, double* data, char * name) {
    int i, j;
    printf("%s[%d, %d]:\n", name, size, size);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            printf("\t%.3f", data[i*size + j]);
        }
        printf("\n");
    }
    printf("\n");
}
