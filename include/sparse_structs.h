#ifndef BBB1CFF6_E369_4CC6_AE40_796811BE778B
#define BBB1CFF6_E369_4CC6_AE40_796811BE778B

#include "mmio.h"

typedef struct {
    int          nrow         ;  // numero de linhas
    int          ncol         ;  // numero de colinhas
    int          nnz          ;  // numero de não zeros           
    int*         row_indices  ;  // indices da linha
    int*         col_indices  ;  // indices da coluna
    double*       values      ;  // coeficientes não nulos
    int is_symmetric          ;  
} crd_t;

typedef struct {
    int nrows;
    int ncols;
    int nnz;
    int *ia;
    int *ja;
    int is_symmetric ;
    double *values;
} csr_t;

crd_t* CRDRead(const char* fname);

void   CRDDestroy(crd_t** crd);

void   CSRDestroy(csr_t** csr);

csr_t* crd2csr(crd_t* data);

void   CSRPrint(csr_t* csr);







#endif /* BBB1CFF6_E369_4CC6_AE40_796811BE778B */
