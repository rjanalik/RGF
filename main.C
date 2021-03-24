#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include "RGF.H"
using namespace std;

#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif

int main(int argc, char *argv[])
{
    int i;
    double data;
    double t0, t1, t2, t3;
    T *b;
    T *x;
    T *invDiag;
    int nrhs;
    RGF<T> *solver;

    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf ("The current date/time is: %s\n",asctime(timeinfo));

    size_t ns = atoi(argv[1]);
    size_t nt = atoi(argv[2]);
    size_t nd = atoi(argv[3]);

    size_t size = ns*nt + nd;

    // load matrix from file
    FILE *F = fopen(argv[4],"r");
   
    size_t fn, fnnz;
    size_t *ia, *ja;
    T *a;
    double val;

    /* read in matrix A, sparse matrix in CSR format */
    fscanf(F,"%zu",&fn);
    fscanf(F,"%zu",&fn);
    fscanf(F,"%zu",&fnnz);

    // allocate memory
    ia = new size_t[fn+1];
    ja = new size_t[fnnz];
    a = new T[fnnz];
  
    for (i = 0; i <= fn; i++){
       fscanf(F,"%zu",&ia[i]);
    }

    for (i = 0; i < ia[fn]; i++){
       fscanf(F,"%zu",&ja[i]);
    }

    for (i = 0; i < ia[fn]; i++){
       fscanf(F,"%lf",&val);
       a[i] = assign_T(val);
    }

    fclose(F);

    nrhs   = 1;
    b      = new T[nrhs*size];
    x      = new T[nrhs*size];
    invDiag= new T[size];

    F = fopen(argv[5],"r");
    for (int i = 0; i < nrhs*size; i++)
    {
       fscanf(F,"%lf",&b[i]);
    }
    fclose(F);
    
    solver = new RGF<T>(ia, ja, a, ns, nt, nd);

    t0 = get_time(0.0);
    solver->factorize();
    t0 = get_time(t0);

    t1 = get_time(0.0);
    solver->solve(x, b, nrhs);
    t1 = get_time(t1);

    t2 = get_time(0.0);
    double logdet = solver->logDet();
    t2 = get_time(t2);

    t3 = get_time(0.0);
    solver->RGFdiag(invDiag);
    t3 = get_time(t3);

    printf("logdet: %f\n", logdet);
    printf("\n");

    printf("factorize time: %lg\n",t0);
    printf("solve time: %lg\n",t1);
    printf("logdet time: %lg\n",t2);
    printf("RGF time: %lg\n",t3);
    printf("\n");

    printf("Residual norm: %e\n", solver->residualNorm(x, b));
    printf("Residual norm normalized: %e\n", solver->residualNormNormalized(x, b));
    printf("\n");

    for (int i = 0; i < nrhs*size; i++)
    {
       printf("%32.24e\n", x[i]);
       //printf("%32.24e\n", invDiag[i]);
    }

    // free memory
    delete[] ia;
    delete[] ja;
    delete[] a;

    delete[] b;
    delete[] x;
    delete solver;
    delete[] invDiag;
    
    return 0;
}

/************************************************************************************************/
