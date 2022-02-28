#ifndef _SPARSE_LINEAR_SOLVER
#define _SPARSE_LINEAR_SOLVER

#include <stdio.h>
#include <stdlib.h>

#include "sparse_structs.h"


double *AllocVector(int n);

/*
 *   Implementa a multiplicação Matriz-Vetor CSR
 */
void CSRMatVec(csr_t* mat, double *x, double *y);

double dot_produt(int n, double *v1, double *v2);

void axpy(int n, double a, double *x, double *y);

void copy(int n, double *vo, double* vd);

void scal(int n, double a, double*v);

void set(int n, double a, double *v);

void print_vector(int n, double*v, const char* name);

void GradientConjugate(csr_t* A, double* b, double *x, double toler, int maxIters, short enable_verbose);

void GradientConjugatePrecond(csr_t* A, double* b, double *x, double toler, int maxIters, short enable_verbose);


#endif


