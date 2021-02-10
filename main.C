#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include "CSR.H"
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
    FILE *F1,*F2;
    int i;
    int size,rank;
    double data;
    double t0;
    TCSR<T> *M;
    T *GR;
    T *b;
    T *invDiag;
    int nrhs;
    RGF<T> *solver;

    if(!rank){
      time_t rawtime;
      struct tm *timeinfo;

      time(&rawtime);
      timeinfo = localtime(&rawtime);
      printf ("The current date/time is: %s\n",asctime(timeinfo));
    }

    int ns = atoi(argv[1]);
    int nt = atoi(argv[2]);
    int nd = atoi(argv[3]);
    

    // load matrix from file
    FILE *F = fopen(argv[4],"r");
   
    int fn, fnnz;
    int *ia, *ja;
    T *a;
    double val;

    /* read in matrix A, sparse matrix in CSR format */
    fscanf(F,"%i",&fn);
    fscanf(F,"%i",&fn);
    fscanf(F,"%i",&fnnz);

    // allocate memory
    ia = new int[fn+1];
    ja = new int[fnnz];
    a = new T[fnnz];
  
    for (i = 0; i <= fn; i++){
       fscanf(F,"%i",&ia[i]);
    }

    for (i = 0; i < ia[fn]; i++){
       fscanf(F,"%i",&ja[i]);
    }

    for (i = 0; i < ia[fn]; i++){
       fscanf(F,"%lf",&val);
       a[i] = assign_T(val);
    }

    fclose(F);

    M      = new TCSR<T>(ia, ja, a, ns, nt, nd);
    GR     = new T[3*(nt-1)*(ns*ns) + (ns*ns)];
    nrhs   = 2;
    b      = new T[nrhs*(ns*nt+nd)];
    invDiag= new T[M->size];

    for (int i = 0; i < nrhs*(ns*nt+nd); i++)
       b[i] = assign_T(i+1);
    
    solver = new RGF<T>(M);

    t0 = get_time(0.0);
    //solver->solve_equation(GR);
    solver->factorize();
    solver->solve(b, nrhs);
    solver->RGFdiag(invDiag);
    t0 = get_time(t0);

    printf("RGF time: %lg\n",t0);

    //// extract diag from GR
    //// for Lisa begin
    //int *GRdiag_ind;
    //T *GRdiag;
    //GRdiag_ind = new int[ns*nt];
    //GRdiag = new T[ns*nt];
    //int ind = 0;
    //for(int i_nt = 0; i_nt < nt; i_nt++)
    //{
    //   for(int i_ns = 0; i_ns < ns-1; i_ns++)
    //   {
    //      int i = i_nt*ns + i_ns;
    //      GRdiag_ind[i] = ind;
    //      ind += ns+1;
    //   }
    //   int i = (i_nt+1)*ns - 1;
    //   GRdiag_ind[i] = ind;
    //   ind++;
    //}
    //for (int i = 0; i < ns*nt; i++)
    //{
    //   GRdiag[i] = GR[GRdiag_ind[i]];
    //}
    //// print diag
    //for (int i = 0; i < ns*nt; i++)
    //{
    //   printf("GRdiag[%d] = %f\n", i, GRdiag[i]);
    //}
    //delete[] GRdiag_ind;
    //delete[] GRdiag;
    //// for Lisa end
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

    delete[] invDiag;
    delete[] GR;
    delete[] b;
    delete M;
    delete solver;
    
    // free memory
    delete[] ia;
    delete[] ja;
    delete[] a;
    
    return 0;
}

/************************************************************************************************/
