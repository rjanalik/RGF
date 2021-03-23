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

/*
Start simulations with RGF NBlock Bmin.dat Bmax.dat M.dat
NBlock: number of blocks of the matrix
Bmin.dat: name of file that contains the first index of each block starting with 0
Bmax.dat: name of file that contains the last index of each block with Bmax[i]=Bmin[i+1]
M.dat: name of file that contains the matrix to work on.
Data stored:
size (matrix size)
n_nonzeros (number of non-zero elements)
fortran index (0 or 1)
index_i index_j real imag (4 columns per matrix entry)
*/

int main(int argc, char *argv[])
{
    int i;
    double data;
    double t0;
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

    nrhs   = 2;
    b      = new T[nrhs*size];
    x      = new T[nrhs*size];
    invDiag= new T[size];

    for (int i = 0; i < nrhs*(ns*nt+nd); i++)
       b[i] = assign_T(i+1);
    
    solver = new RGF<T>(ia, ja, a, ns, nt, nd);

    t0 = get_time(0.0);

    solver->factorize();
    solver->solve(x, b, nrhs);

    printf("logdet: %f\n", solver->logDet());

    solver->RGFdiag(invDiag);
    t0 = get_time(t0);

    printf("RGF time: %lg\n",t0);

    printf("Residual norm: %e\n", solver->residualNorm(x, b));
    printf("Residual norm normalized: %e\n", solver->residualNormNormalized(x, b));

    for (int i = 0; i < nrhs*(ns*nt+nd); i++)
    {
       printf("x[%d] = %f\n", i, b[i]);
    }
    printf("\n");
    //for (int i = 0; i < M->size; i++)
    for (int i = 0; i < M->size; i++)
    {
       printf("invDiag[%d] = %f\n", i, invDiag[i]);
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
