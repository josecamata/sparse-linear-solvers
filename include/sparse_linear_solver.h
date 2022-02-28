#ifndef _SPARSE_LINEAR_SOLVER
#define _SPARSE_LINEAR_SOLVER

#include <stdio.h>
#include <stdlib.h>

#include "sparse_structs.h"


double *AllocVector(int n);

/**
 * @brief Produto matriz vetor para matrizes esparsas
 *        no formato CSR
 *        y = A*x 
 * 
 * @param mat matrix CSR
 * @param x   vetor que multiplica mat
 * @param y   vetro com o resultante
 */
void CSRMatVec(csr_t* mat, double *x, double *y);

/**
 * @brief produto interno entre dois vetores
 *        dot_produt = v1*v2'
 * 
 * @param n  tamanho do vetores
 * @param v1 vetor 1
 * @param v2 vetor 2
 * @return double 
 */
double dot_produt(int n, double *v1, double *v2);

/**
 * @brief Operação y = y + ax
 * 
 * @param n tamanho dos vetores x e y
 * @param a escalar 
 * @param x vetor
 * @param y vetor resultante
 */
void axpy(int n, double a, double *x, double *y);

/**
 * @brief copaia o conteudo de vo para vd;
 * 
 * @param n  tamanha dos vetores
 * @param vo vetor de origem
 * @param vd vetor de destino
 */
void copy(int n, double *vo, double* vd);

/**
 * @brief Multiplica um vetor por um escalar
 *        v = a*v
 * @param n Tamanho do Vetor
 * @param v Vetor
 * @param a Escalar
 */
void scal(int n, double*v, double a);

/**
 * @brief  Atribui uma valor constante as entradas de um vetor
 *           v[i] = a, para i=0,1,...,n-1
 * @param n Tamanho do vetor
 * @param v Vetor
 * @param a Escalar
 */
void set(int n, double *v, double a);

/**
 * @brief Imprime conteudo do vetor na tela
 * 
 * @param n Tamanho do vetor c
 * @param v Vetor de reais
 * @param name Identidicação do vetor
 */
void print_vector(int n, double*v, const char* name);

/**
 * @brief Solucionador Gradiente Conjugado empregado na
 *        solução de sistemas simetricos
 * 
 * @param A Matriz CSR
 * @param b Vetor dos termos independente
 * @param x Vetor solução
 * @param toler precisão empregada
 * @param maxIters numero máximo de iterações
 * @param enable_verbose habita impressão do residuo por iteração
 */
void GradientConjugate(csr_t* A, double* b, double *x, double toler, int maxIters, short enable_verbose);

/**
 * @brief Solucionador Gradiente Conjugad empregado na
 *        solução de sistemas simetricos.
 *        Usa precondicionador diagonal.
 * 
 * @param A Matriz CSR
 * @param b Vetor dos termos independente
 * @param x Vetor solução
 * @param toler precisão empregada
 * @param maxIters numero máximo de iterações
 * @param enable_verbose habita impressão do residuo por iteração
 */
void GradientConjugatePrecond(csr_t* A, double* b, double *x, double toler, int maxIters, short enable_verbose);


#endif


