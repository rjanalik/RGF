#ifndef __CSR
#define __CSR

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "Types.H"

template <class T>
class TCSR{
public:
    
    TCSR(int,int,int);
    TCSR(char*);
    TCSR(int*,int*,T*,int,int);
    ~TCSR();

    void get_row_edge();

    int pos,size,size_tot,n_nonzeros,type,findx,first_row;
    T *nnz;
    int *index_i,*index_j,*edge_i,*diag_pos;

    int NBlock;
    int *Bmin,*Bmax;
    
private:
};

template <class T>
TCSR<T>::TCSR(int N, int n_nnz, int fortran_index)
{
    nnz        = new T[max(n_nnz,1)];
    index_i    = new int[max(N,1)];
    index_j    = new int[max(n_nnz,1)];
    edge_i     = new int[N+1];
    diag_pos   = new int[max(N,1)];

    size       = N;
    size_tot   = N;
    type       = 0;
    n_nonzeros = n_nnz;
    findx      = fortran_index;
    first_row  = 0;                                                    
}

/************************************************************************************************/

template <>
inline TCSR<CPX>::TCSR(char *filename)
{

    int i,ind_i;
    double r_nnz,i_nnz;
    FILE *F = fopen(filename,"r");
   
    fscanf(F,"%i",&size);
    fscanf(F,"%i",&n_nonzeros);
    fscanf(F,"%i",&findx);
    
    nnz        = new CPX[n_nonzeros];
    index_i    = new int[size];                                           
    index_j    = new int[n_nonzeros];                                       
    edge_i     = new int[size+1];                                         
    diag_pos   = new int[size];                                           

    size_tot   = size;                                                    
    type       = 0;                                        
    first_row  = 0;                                                    

    for(i=0;i<size;i++) index_i[i] = 0;
   
    for(i=0;i<n_nonzeros;i++){

        fscanf(F,"%i",&ind_i);
	fscanf(F,"%i",index_j+i);
	fscanf(F,"%lg",&r_nnz);
	fscanf(F,"%lg",&i_nnz);
  
	nnz[i] = CPX(r_nnz,i_nnz);

	if(ind_i==index_j[i]){
	    diag_pos[ind_i-findx] = i;
	}
	index_i[ind_i-findx]++;
	
    }
    get_row_edge();
    fclose(F);

}

// RJ: load CSR, make dense blocks
template <class T>
inline TCSR<T>::TCSR(int *ia, int *ja, T *a, int ns, int nt)
{

    int i,ind_i;
    double r_nnz,i_nnz;
   
    size = ns*nt;
    n_nonzeros = 2*(nt-1)*ns*ns + ns*ns;
    findx = 0; // 0-C index, 1-fortran index
    
    nnz        = new T[n_nonzeros];
    index_i    = new int[size];                                           
    index_j    = new int[n_nonzeros];                                       
    edge_i     = new int[size+1];                                         
    diag_pos   = new int[size];                                           

    size_tot   = size;                                                    
    type       = 0;                                        
    first_row  = 0;                                                    

    NBlock = nt;
    Bmin   = new int[NBlock];
    Bmax   = new int[NBlock];
    for (int i = 0; i < nt; i++)
    {
       Bmin[i] = i*ns;
       Bmax[i] = (i+1)*ns;
    }

    for(i=0;i<size;i++) index_i[i] = 0;

    // create a sparse matrix with known structure
    // initialize with 0.0
    int i_j = 0;
    for(int i_nt = 0; i_nt < nt-1; i_nt++)
    {
       for(int i_ns = 0; i_ns < ns; i_ns++)
       {
          int i_r = i_nt*ns + i_ns;
          int i_c_f = i_nt*ns;
          int i_c_l = (i_nt+2)*ns;
          edge_i[i_r] = i_j;
          for(int i_c = i_c_f; i_c < i_c_l; i_c++)
          {
             index_j[i_j++] = i_c;
          }
       }
    }
    for(int i_ns = 0; i_ns < ns; i_ns++)
    {
       int i_r = (nt-1)*ns + i_ns;
       int i_c_f = (nt-1)*ns;
       int i_c_l = size;
       edge_i[i_r] = i_j;
       for(int i_c = i_c_f; i_c < i_c_l; i_c++)
       {
          index_j[i_j++] = i_c;
       }
    }
    edge_i[size] = i_j;
    for(i = 0; i < n_nonzeros; i++)
    {
       nnz[i] = T(0);
    }
    for(int r = 0; r < size; r++)
    {
       int i = edge_i[r];
       while(index_j[i] < r)
          i++;
       diag_pos[r] = i;
    }
    for(i = 0; i < size; i++)
    {
       index_i[i] = edge_i[i+1] - edge_i[i];
    }
    //get_row_edge();

    // add the matrices (overwrite, we know there are zeros)
    for(int i_r = 0; i_r < size; i_r++)
    {
       int i = edge_i[i_r];
       for(int fi_c = ia[i_r]; fi_c < ia[i_r+1]; fi_c++)
       {
          int c = ja[fi_c];
          while(index_j[i] < c)
             i++;
          nnz[i]   = a[fi_c];
       }
    }

}

/***********************************************************************************************/

template <>
inline TCSR<double>::TCSR(char *filename)
{

    int i,ind_i;
    FILE *F = fopen(filename,"r");
   
    fscanf(F,"%i",&size);
    fscanf(F,"%i",&n_nonzeros);
    fscanf(F,"%i",&findx);
    
    nnz        = new double[n_nonzeros];
    index_i    = new int[size];                                           
    index_j    = new int[n_nonzeros];                                       
    edge_i     = new int[size+1];                                         
    diag_pos   = new int[size];                                           

    size_tot   = size;                                                    
    type       = 0;                                        
    first_row  = 0;                                                    

    for(i=0;i<size;i++) index_i[i] = 0;
   
    for(i=0;i<n_nonzeros;i++){
        fscanf(F,"%i",&ind_i);
	fscanf(F,"%i",index_j+i);
	fscanf(F,"%lg",nnz+i);
	if(ind_i==index_j[i]){
	    diag_pos[ind_i-findx] = i;
	}
	index_i[ind_i-findx]++;
	
    }
    get_row_edge();
    fclose(F);

}

/***********************************************************************************************/

template <class T>
TCSR<T>::~TCSR()
{
    delete[] nnz;
    delete[] index_i;
    delete[] index_j;
    delete[] edge_i;
    delete[] diag_pos;

    delete[] Bmin;
    delete[] Bmax;
}

/************************************************************************************************/

template <class T>
void TCSR<T>::get_row_edge()
{
    int i;

    edge_i[0] = findx;
  
    for(i=0;i<size;i++){
        edge_i[i+1] = edge_i[i]+index_i[i];
    }
    
}

/************************************************************************************************/

#endif
