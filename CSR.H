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
    TCSR(size_t*,size_t*,T*,size_t,size_t,size_t);
    ~TCSR();

    void get_row_edge();

    size_t pos,size,size_tot,n_nonzeros,type,findx,first_row;
    T *nnz;
    size_t *index_i,*index_j,*edge_i,*diag_pos;

    size_t NBlock;
    size_t *Bmin,*Bmax;
    size_t ns, nt, nd;
    
private:
};

template <class T>
TCSR<T>::TCSR(int N, int n_nnz, int fortran_index)
{
    nnz        = new T[max(n_nnz,1)];
    index_i    = new size_t[max(N,1)];
    index_j    = new size_t[max(n_nnz,1)];
    edge_i     = new size_t[N+1];
    diag_pos   = new size_t[max(N,1)];

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
    index_i    = new size_t[size];
    index_j    = new size_t[n_nonzeros];
    edge_i     = new size_t[size+1];
    diag_pos   = new size_t[size];

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
inline TCSR<T>::TCSR(size_t *ia, size_t *ja, T *a, size_t _ns, size_t _nt, size_t _nd)
{

    size_t i,ind_i;
    double r_nnz,i_nnz;

    ns = _ns;
    nt = _nt;
    nd = _nd;
   
    size = ns*nt + nd;
    n_nonzeros = 2*(nt-1)*ns*ns + ns*ns + size*nd;
    findx = 0; // 0-C index, 1-fortran index
    
    nnz        = new T[n_nonzeros];
    index_i    = new size_t[size];
    index_j    = new size_t[n_nonzeros];
    edge_i     = new size_t[size+1];
    diag_pos   = new size_t[size];

    size_tot   = size;                                                    
    type       = 0;                                        
    first_row  = 0;                                                    

    if (nd > 0)
       NBlock = nt + 1;
    else
       NBlock = nt;
    Bmin   = new size_t[NBlock];
    Bmax   = new size_t[NBlock];
    for (size_t i = 0; i < nt; i++)
    {
       Bmin[i] = i*ns;
       Bmax[i] = (i+1)*ns;
    }
    if (nd > 0)
    {
       Bmin[nt] = nt*ns;
       Bmax[nt] = nt*ns + nd;
    }

    for(i=0;i<size;i++) index_i[i] = 0;

    // create a sparse matrix with known structure
    // initialize with 0.0
    size_t i_j = 0;
    for(size_t i_nt = 0; i_nt < nt-1; i_nt++)
    {
       for(size_t i_ns = 0; i_ns < ns; i_ns++)
       {
          size_t i_r = i_nt*ns + i_ns;
          size_t i_c_f = i_nt*ns;
          size_t i_c_l = (i_nt+2)*ns;
          edge_i[i_r] = i_j;
          for(size_t i_c = i_c_f; i_c < i_c_l; i_c++)
          {
             index_j[i_j++] = i_c;
          }
          if (nd > 0)
          {
             for(size_t i_c = Bmin[nt]; i_c < Bmax[nt]; i_c++)
             {
                index_j[i_j++] = i_c;
             }
          }
       }
    }
    for(size_t i_ns = 0; i_ns < ns; i_ns++)
    {
       size_t i_r = (nt-1)*ns + i_ns;
       size_t i_c_f = (nt-1)*ns;
       size_t i_c_l = size;
       edge_i[i_r] = i_j;
       for(size_t i_c = i_c_f; i_c < i_c_l; i_c++)
       {
          index_j[i_j++] = i_c;
       }
    }
    for(size_t i_nd = 0; i_nd < nd; i_nd++)
    {
       size_t i_r = Bmin[nt] + i_nd;
       edge_i[i_r] = i_j;
       for(size_t i_c = Bmin[nt]; i_c < Bmax[nt]; i_c++)
       {
          index_j[i_j++] = i_c;
       }
    }
    edge_i[size] = i_j;
    for(i = 0; i < n_nonzeros; i++)
    {
       nnz[i] = T(0);
    }
    for(size_t r = 0; r < size; r++)
    {
       size_t i = edge_i[r];
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
    for(size_t i_r = 0; i_r < size; i_r++)
    {
       size_t i = edge_i[i_r];
       for(size_t fi_c = ia[i_r]; fi_c < ia[i_r+1]; fi_c++)
       {
          size_t c = ja[fi_c];
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
    index_i    = new size_t[size];
    index_j    = new size_t[n_nonzeros];
    edge_i     = new size_t[size+1];
    diag_pos   = new size_t[size];

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
