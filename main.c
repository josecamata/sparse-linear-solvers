/*
    Programa para demonstrar a método
    do gradiente conjugado

    Autor: Jose Camata
    E-mail: camata@ice.ufjf.br 

*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "sparse_structs.h"
#include "sparse_linear_solver.h"


/* Mostra a ajuda */
void show_help(char *name) {
    fprintf(stderr, "\
    [uso] %s <opcoes>\n\
            -h                   mostra essa tela e sai.\n\
            -m [ARQUIVO .mtx]    arquivo Matrix Market  .\n\
            -v                   Imprime saida completa .\n", name) ;
    exit(-1) ;
}

// Define um enumerador para listar as opções de impressão na tela
// VERBOSE_OFF: imprime informação basica
// VERBOSE_ON: Imprime o residuo por iteração
typedef enum {VERBOSE_OFF=0, VERBOSE_ON=1} verbose_mode_t;

int main(int argc, char *argv[])
{

    int opt;

     /* Variáveis que receberão os argumentos
     * das opções. */
    char *file=NULL;
    short flag=0;
    verbose_mode_t  verbose_mode = VERBOSE_OFF;

    /* Chama ajuda. */
    if ( argc < 2 ) show_help(argv[0]) ;

    /* getopt() retorna o caractere de uma opção a cada
     * iteração e -1 para marcar o fim do processo. */
    while( (opt = getopt(argc, argv, "hvm:")) > 0 ) {
        switch ( opt ) {
            case 'h': /* help */
                show_help(argv[0]) ;
                break ;
            case 'm': /* opção -m */
                file = optarg ; // optarg é uma variavel que aponta para o nome do arquivo
                flag =1;
                break ;
            case 'v': /* opção -v */
                verbose_mode = VERBOSE_ON;
                break ;
            default:
                fprintf(stderr, "Opcao invalida ou faltando argumento: `%c'\n", opt) ;
                return -1 ;
        }
    }
    if(!flag)
    {
        fprintf(stderr, "Arquivo de dados não especificado\n") ;
        show_help(argv[0]);    
        return -1;
    }

    printf("============================================\n");
    printf("|                 DCC089                   |\n");
    printf("|           Topicos Especias em            |\n");
    printf("|                  em                      |\n");
    printf("|          Computação Cientifica I         |\n");
    printf("|------------------------------------------|\n");
    printf("|  Métodos dos Gradientes Conjugados       |\n");
    printf("============================================\n");


    // Cria a estrutura para armazenamento da matriz 
    // no formanto coordenadas;
    crd_t *crd;
    crd        = CRDRead(file);

    // Converte o formato coordenada para CSR
    csr_t *A   = crd2csr(crd);
    
    // Apenas imprime a matriz se ela for pequena.
    if(verbose_mode && A->nrows < 10)
        CSRPrint(A);

    // Verifica se a matriz é quadrada.
    if(A->ncols == A->nrows)
    {
        int n     = A->nrows;

        printf("Numero que incognitas: %d\n",n);
        double *u = AllocVector(n);
        double *x = AllocVector(n);
        double *b = AllocVector(n);
        double *e = AllocVector(n);

        // Solução Exata
        // u = 1.0;
        set(n,1.0,u);
        
        // Solução aproximada inicial
        // x = 0.0
        set(n,0.0,x);

        // Lado direito é obtido multiplicado A pela solução exata u
        // b = A*u;
        CSRMatVec(A,u,b);

        // Resolve o sistema Ax = b usando o método do Gradiente Conjugado
        GradientConjugate(A,b,x,1.0E-6,n,verbose_mode);

        // Verifica a solução x tem que ser aproximadamente igual a u
        copy(n,u,e);
        axpy(n,-1.0,x,e);
        printf("|U_exato - U_aprox| = %e\n" , sqrt(dot_produt(n,e,e)));

   
        // Reescreve solução aproximada inicial
        // x = 0.0
        set(n,0.0,x);

        // Resolve o sistema Ax = b usando o método do 
        // Gradiente Conjugado com precondicionador Jacobi
        GradientConjugatePrecond(A,b,x,1.0E-6,n,verbose_mode);

        // Verifica a solução
        copy(n,u,e);
        axpy(n,-1.0,x,e);
        printf("|U_exato - U_aprox| = %e\n" , sqrt(dot_produt(n,e,e)));
        
        free(x);
        free(b);
        free(u);
    }
    CRDDestroy(&crd);
    CSRDestroy(&A);
    return 0;
}