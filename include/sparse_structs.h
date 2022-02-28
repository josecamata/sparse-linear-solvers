#ifndef BBB1CFF6_E369_4CC6_AE40_796811BE778B
#define BBB1CFF6_E369_4CC6_AE40_796811BE778B

#include "mmio.h"

/**
 * @brief Estrutura de dados para armazenar dados do formato
 *        Matrix Market
 */
typedef struct {
    int          nrow         ;  // numero de linhas 
    int          ncol         ;  // numero de colinhas
    int          nnz          ;  // numero de não zeros           
    int*         row_indices  ;  // indices da linha
    int*         col_indices  ;  // indices da coluna
    double*       values      ;  // coeficientes não nulos
    int is_symmetric          ;  
} crd_t;

/**
 * @brief Estrutura de dados para armazenar matriz esparsa no 
 *        formato CSR
 */
typedef struct {
    int nrows;
    int ncols;
    int nnz;
    int *ia;
    int *ja;
    int is_symmetric ;
    double *values;
} csr_t;

/**
 * @brief  Le arquivo no formanto matrix market
 * 
 * @param fname no arquivo .mtx
 * @return crd_t* :Ponteiro para a estrutura de dados 
 */
crd_t* CRDRead(const char* fname);

/**
 * @brief Desaloca memoeria da estrutura
 * 
 * @param crd 
 */
void   CRDDestroy(crd_t** crd);

void   CSRDestroy(csr_t** csr);

/**
 * @brief Converte CRD para CSR
 * 
 * @param data  Estrutura CRD
 * @return csr_t*:
 */
csr_t* crd2csr(crd_t* data);

void   CSRPrint(csr_t* csr);


#endif /* BBB1CFF6_E369_4CC6_AE40_796811BE778B */
