#include <mpi.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include "CSR.H"
#include "RGF.H"
using namespace std;

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
    int NBlock;
    int Bsize;
    int *Bmin,*Bmax;
    double data;
    double t0;
    TCSR<CPX> *M;
    CPX *GR;
    RGF *solver;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(!rank){
      time_t rawtime;
      struct tm *timeinfo;

      time(&rawtime);
      timeinfo = localtime(&rawtime);
      printf ("The current date/time is: %s\n",asctime(timeinfo));
    }

    NBlock = atoi(argv[1]);

    F1     = fopen(argv[2],"r");
    F2     = fopen(argv[3],"r");

    Bmin   = new int[NBlock];
    Bmax   = new int[NBlock];
    Bsize  = 0;
    
    for(i=0;i<NBlock;i++){
      
        fscanf(F1,"%lg",&data);
	Bmin[i] = Round(data);
	
	fscanf(F2,"%lg",&data);
	Bmax[i] = Round(data);
	
	if((Bmax[i]-Bmin[i])>Bsize){
	    Bsize = Bmax[i]-Bmin[i];
	}
    }

    fclose(F1);
    fclose(F2);

    M      = new TCSR<CPX>(argv[4]);
    GR     = new CPX[Bmax[NBlock-1]*Bsize];
    
    solver = new RGF(M);

    t0 = get_time(0.0);
    solver->solve_equation(GR,Bmin,Bmax,NBlock);
    t0 = get_time(t0);

    printf("RGF time: %lg\n",t0);
    
    delete[] Bmin;
    delete[] Bmax;
    delete[] GR;
    delete M;
    delete solver;
    
    MPI_Finalize();
    
    return 0;
}

/************************************************************************************************/
