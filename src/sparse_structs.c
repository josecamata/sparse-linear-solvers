
#include <stdlib.h>
#include <stdio.h>

#include "sparse_structs.h"


crd_t* CRDRead(const char* fname)  // .mtx file
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    crd_t *mm_matrix;

    int         nrow ;  // numero de linhas
    int         ncol ;  // numero de colinhas
    int          nnz ;  // numero de não zeros           

    if ((f = fopen(fname, "r")) == NULL) {
        printf("Arquivo não encontrado: %s\n", fname);
        exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f,&nrow,&ncol,&nnz)) !=0)
        exit(1);

    mm_matrix= (crd_t*) malloc (sizeof(crd_t));
    /* reseve memory for matrices */
    mm_matrix->nrow = nrow;
    mm_matrix->ncol = ncol;
    mm_matrix->nnz  = nnz;
    mm_matrix->row_indices = (int *) malloc(nnz * sizeof(int));
    mm_matrix->col_indices = (int *) malloc(nnz * sizeof(int));
    mm_matrix->values      = (double *) malloc(nnz * sizeof(double));
    mm_matrix->is_symmetric =  mm_is_symmetric(matcode);

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (int i=0; i<nnz; i++)
    {
        fscanf(f, "%d %d %lg\n", &mm_matrix->row_indices[i], &mm_matrix->col_indices[i], &mm_matrix->values[i]);
        mm_matrix->row_indices[i]--;  /* adjust from 1-based to 0-based */
        mm_matrix->col_indices[i]--;
    }

    if (f !=stdin) fclose(f);

    return mm_matrix;

}


void CRDDestroy(crd_t **crd)
{
    free((*crd)->col_indices);
    free((*crd)->row_indices);
    free((*crd)->values);
    free((*crd));
}


void CSRDestroy(csr_t **crd)
{
    free((*crd)->ia);
    free((*crd)->ja);
    free((*crd)->values);
    free((*crd));
}

/*
 *
 *
 * 
 * */
csr_t* crd2csr(crd_t* data)
{

    csr_t* sparse = (csr_t*) malloc (sizeof(csr_t));

    sparse->nrows = data->nrow;
    sparse->ncols = data->ncol;  
    sparse->is_symmetric = data->is_symmetric;

    sparse->ia       = (int*) malloc((sparse->nrows+1)*sizeof(int));

    int *nnz_per_row = (int*) malloc((sparse->nrows)*sizeof(int));

    for(int i =0; i <= sparse->nrows; ++i)
        sparse->ia[i] = 0;
    // zerar nnz_per_row;
    for(int i =0; i < sparse->nrows; i++)
        nnz_per_row[i] = 0;
    
    // Conta o numero de não zeros por linha
    int nnz_total = 0;
    for(int idx =0; idx < data->nnz; ++idx)
    {
        int i = data->row_indices[idx];
        int j = data->col_indices[idx];

        if(!data->is_symmetric) {
            nnz_per_row[i]++;
            nnz_total++;
        } else
        {
                
            nnz_per_row[j]++;
            nnz_total++;
        }
        
    }

    sparse->ja     = (int*)    malloc((nnz_total)*sizeof(int));
    sparse->values = (double*) malloc((nnz_total)*sizeof(double));

    // Calcula o offset
    sparse->ia[0] = 0;
    for(int i =1; i <= sparse->nrows; i++)
        sparse->ia[i] = sparse->ia[i-1]+nnz_per_row[i-1];

    // zerar nnz_per_row;
    for(int i =0; i < sparse->nrows; i++)
        nnz_per_row[i] = 0;


    for(int idx =0; idx < data->nnz; ++idx)
    {
        int i          = data->row_indices[idx];
        int j          = data->col_indices[idx];
        double val     = data->values[idx];

        if(!data->is_symmetric)
        {
            int    offset  = sparse->ia[i] + nnz_per_row[i];
            sparse->ja[offset]        = j;
            sparse->values[offset]    = val;
            nnz_per_row[i]++;
        }
        else
        {
            int    offset  = sparse->ia[j] + nnz_per_row[j];
            sparse->ja[offset]        = i;
            sparse->values[offset]    = val;
            nnz_per_row[j]++;            
        }
        
    }

    free(nnz_per_row);
    return sparse;

}

void CSRPrint(csr_t *c)
{
    if(c->is_symmetric)
        printf("Matrix (CSR Symmetric) \n");
    else
        printf("Matrix (CSR General) \n");
    for(int r = 0; r < c->nrows; ++r)
    {
        printf("row %d: ", r);
        int jstart = c->ia[r];
        int jend   = c->ia[r+1];
        for(int j=jstart; j < jend; ++j)
            printf("(%d,%f) " , c->ja[j], c->values[j]);
        printf("\n");
    }
}
