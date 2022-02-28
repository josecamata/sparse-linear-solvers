#include "config.h"
#include "sparse_linear_solver.h"
#include <math.h>

#if defined(HAVE_CBLAS)
#include "cblas.h"
#endif

double *AllocVector(int n)
{
    return (double*) malloc (n*sizeof(double));
}

void CSRMatVec(csr_t* mat, double *x, double *y)
{
    if(mat->is_symmetric)
    {
        for(int r = 0; r < mat->nrows; ++r)
        {
            int jstart = mat->ia[r];
            int jend   = mat->ia[r+1];
            y[r] = 0.0;
            // Multiplica por (D+U)
            for(int j=jstart; j < jend; ++j) {
                int c = mat->ja[j];
                y[r] += mat->values[j]*x[c];
            }
        }
        for(int r = 0; r < mat->nrows; ++r)
        {
            int jstart = mat->ia[r];
            int jend   = mat->ia[r+1];
            // Multipluca por L
            for(int j=jstart+1; j < jend; ++j) {
                int c = mat->ja[j];
                y[c] += mat->values[j]*x[r];
            }
        } 
    }
    else
    {
        for(int r = 0; r < mat->nrows; ++r)
        {
            int jstart = mat->ia[r];
            int jend   = mat->ia[r+1];
            double s = 0.0;
            for(int j=jstart; j < jend; ++j)
                s += mat->values[j]*x[mat->ja[j]];
            y[r] = s; 
        }
    }
}

double dot_produt(int n, double *v1, double *v2)
{
#if defined(HAVE_CBLAS)
    return cblas_ddot(n,v1,1,v2,1);
#else
    double r = 0.0;
    for(int i = 0; i < n; i++)
        r += v1[i]*v2[i];
    return r;
#endif
}


void axpy(int n, double a, double *x, double *y)
{
#if defined(HAVE_CBLAS)
    cblas_daxpy(n,a,x,1,y,1);
#else
    for(int i = 0; i < n; i++)
    {
        y[i]+= a*x[i];
    }
#endif
}

void copy(int n, double *vo, double* vd)
{
#if defined(HAVE_CBLAS)
    cblas_dcopy(n,vo,1,vd,1);
#else
    for(int i = 0; i < n; i++)
        vd[i] = vo[i];
#endif
}

void scal(int n, double a, double *v)
{
#if defined(HAVE_CBLAS)
    cblas_dscal(n,a,v,1);
#else
    for(int i = 0; i < n; i++)
         v[i] = a*v[i];
#endif
}

void set(int n, double a, double *v)
{
    for(int i = 0; i < n; i++)
        v[i] = a;

}

void GradientConjugate(csr_t* A, double* b, double *x, double toler, int maxIters,short enable_verbose)
{
    int n     = A->nrows;
    // working vector
    double *p = AllocVector(n);
    double *v = AllocVector(n);
    double *r = AllocVector(n);
    double rdot_old, rdot_new, tmp;
    double error, alpha, beta;
    int k = 0;

    printf("Gradiente Conjugado\n");
    // r = b
    copy(n,b,r);
    // v = A*x
    CSRMatVec(A,x,v);
    // r = r - v
    axpy(n,-1.0,v,r);
    // rsold = |r|
    rdot_old = dot_produt(n,r,r);

    // p = r
    copy(n,r,p);

    while((rdot_old > toler*toler) && (k < maxIters) )
    {
        k = k + 1;

        // v = A*p
        CSRMatVec(A,p,v);

        tmp = dot_produt(n,v,p);
        // alpha = rsold/(v'*p);
        alpha = rdot_old/tmp;

        //x     = x + alpha*p;
        axpy(n,alpha,p,x);

        // r     = r - alpha*v;
        axpy(n,-alpha,v,r);

        //rsnew = r'*r;
        rdot_new = dot_produt(n,r,r);
        
        beta = rdot_new/rdot_old;

        //p    = r + beta*p

        scal(n,beta,p);
        axpy(n,1.0,r,p);

        rdot_old = rdot_new;

        if(enable_verbose)
            printf("iter = %d -- |r| = %f\n", k, sqrt(rdot_old));

    }

   
    printf("Iterações: %d \n", k);
    printf("Residuo final: %f\n", rdot_old);

    free(p);
    free(v);
    free(r);
    
}

//
// Aplica o precondicionar r = inv(M)*v
void apply_precond_jacobi(int n, double* invM, double* v, double *r)
{
    for(int i=0; i < n; ++i)
    {
        r[i] = invM[i]*v[i];
    }
}

// Extrai a matriz precondicionadora.
// Observe que a função ja calcula a inversa da matriz diagonal.
void GetPrecondJacobi(csr_t* mat, double *invM)
{
    for(int r = 0; r < mat->nrows; ++r)
    {
        int jstart = mat->ia[r];
        int jend   = mat->ia[r+1];
        for(int j=jstart; j < jend; ++j)
            if(mat->ja[j] == r) {
                invM[r] = (mat->values[j]!=0.0)?1.0/mat->values[j]:1.0;
                continue;
            }
    }
}


void GradientConjugatePrecond(csr_t* A, double* b, double *x, double toler, int maxIters, short enable_verbose)
{
    int n     = A->nrows;
    // working vectors
    double *p = AllocVector(n);
    double *v = AllocVector(n);
    double *r = AllocVector(n);
    double *z = AllocVector(n);
    double *invM = AllocVector(n);

    double rdot_old, rdot_new, tmp;
    double rdot, alpha, beta;
    int k = 0;

    printf("Gradiente Conjugado Precondicionado\n");
    GetPrecondJacobi(A,invM);

    // r = b
    copy(n,b,r);

    // v = A*x
    CSRMatVec(A,x,v);

    // r = r - v
    axpy(n,-1.0,v,r);

    // z = inv(M)*r
    apply_precond_jacobi(n,invM,r,z);

    // rsold = |r|
    rdot_old = dot_produt(n,r,z);
    // p = z
    copy(n,z,p);

    rdot = dot_produt(n,r,r);
    //
    while((rdot > toler*toler) && (k < maxIters) )
    {
        // 
        k = k + 1;
        // v = A*p
        CSRMatVec(A,p,v);

        tmp = dot_produt(n,v,p);
        // alpha = rsold/(v'*p);
        alpha = rdot_old/tmp;

        //x     = x + alpha*p;
        axpy(n,alpha,p,x);

        // r     = r - alpha*v;
        axpy(n,-alpha,v,r);

        rdot = dot_produt(n,r,r);

        apply_precond_jacobi(n,invM,r,z);

        //rsnew = r'*r;
        rdot_new = dot_produt(n,r,z);
  
        beta = rdot_new/rdot_old;

        //p    = r + beta*p
        scal(n,beta,p);
        axpy(n,1.0,z,p);

        rdot_old = rdot_new;

        if(enable_verbose)
            printf("iter = %d -- |r| = %f\n", k, sqrt(rdot));


    }

   
    printf("Iterações: %d \n", k);
    printf("Residuo final: %f\n", sqrt(rdot));

    free(p);
    free(v);
    free(r);
    free(invM);
    free(z);
    
}

void print_vector(int n, double*v, const char* name)
{
    printf("Vector: %s\n", name);
    for(int i = 0; i < n; i++)
    {
        printf("%d: %f\n", i, v[i]);
    }
}


