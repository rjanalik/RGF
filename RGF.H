#ifndef __RGF
#define __RGF

#include <string.h>

#include "cublas_v2.h"
#include "CSR.H"
#include "Types.H"
#include "Utilities.H"
#include "CWC_utility.H"

template <class T>
class RGF{
	
public:

    RGF(TCSR<T>*);

    ~RGF();
    void solve_equation(T*);
    void factorize();
    void solve(T*,int);
    void solve(T*,T*,int);
    void RGFdiag(T*);
				
private:

    TCSR<T>* matrix;
    int findx;
    int b_size;
    bool MF_dev_allocated = false;

    cudaStream_t stream_c,stream_m;
    void *cublas_handle;
    magma_queue_t magma_queue;
    T *M_dev;
    T *M1_dev;
    T *M2_dev;
    T *M3_dev;
    T *MF_dev;
    T *rhs_dev;
    T *nnz_dev;
    int *edge_i_dev;
    int *index_j_dev;

    inline int mf_block_index(int, int);
    inline int mf_block_lda(int, int);
    inline int mf_dense_block_index(int);
    inline int mf_dense_block_lda(int);

    void extract_diag(int*,int*,T*,T*,int*,int*,int,cudaStream_t);
    void extract_not_diag(int*,int*,T*,T*,int*,int*,int,int,cudaStream_t);
    void FirstStage(int*,int*,int);
    void SecondStage(T*,int*,int*,int);
    void FirstStageFactor(int*,int*,int);
    void SecondStageSolve(int*,int*,int,int);
    void ThirdStageRGF(T*,T*,int*,int*,int);
    void create_blocks(int*,int*,int);
    void write_matrix(const char*,T*,int,int);

   T f_one();
   T f_zero();
};

/************************************************************************************************/

template <class T>
RGF<T>::RGF(TCSR<T>* mat)
{
    matrix = mat;
    findx  = 0;
}

/************************************************************************************************/

template <class T>
RGF<T>::~RGF()
{
   if (MF_dev_allocated)
   {
      deallocate_data_on_dev(MF_dev,matrix->n_nonzeros*sizeof(T));
   }
}

/************************************************************************************************/

template <class T>
inline int RGF<T>::mf_block_index(int r, int c)
{
   return matrix->diag_pos[c*b_size] + (r-c)*b_size;
}

/************************************************************************************************/

template <class T>
inline int RGF<T>::mf_block_lda(int r, int c)
{
   return matrix->index_i[c*b_size];
}

/************************************************************************************************/

template <class T>
inline int RGF<T>::mf_dense_block_index(int i)
{
   if (i < matrix->nt-1)
      return matrix->diag_pos[i*b_size] + 2*b_size;
   else if (i == matrix->nt-1)
      return matrix->diag_pos[i*b_size] + b_size;
   else
      return matrix->diag_pos[i*b_size];
}

/************************************************************************************************/

template <class T>
inline int RGF<T>::mf_dense_block_lda(int i)
{
   return mf_block_lda(i, i);;
}

/************************************************************************************************/

template <class T>
void RGF<T>::extract_diag(int *edge_i,int *index_j,T *nnz,T *D,int *Bmin,int *Bmax,\
		       int index,cudaStream_t stream)
{
    int NR;
    int imin,imax;

    NR   = Bmax[index]-Bmin[index];
    imin = Bmin[index];
    imax = Bmax[index];

    init_var_on_dev(D,NR*NR,stream);

    extract_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,matrix->first_row,matrix->findx,\
			  stream);
}

/************************************************************************************************/

template <class T>
void RGF<T>::extract_not_diag(int *edge_i,int *index_j,T *nnz,T *D,int *Bmin,int *Bmax,\
			   int index,int side,cudaStream_t stream)
{
    int NR,NC;
    int imin,imax;
    int jmin;

    NR   = Bmax[index]-Bmin[index];
    NC   = Bmax[index+side]-Bmin[index+side];
    imin = Bmin[index];
    imax = Bmax[index];
    jmin = Bmin[index+side];

    init_var_on_dev(D,NR*NC,stream);

    extract_not_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,jmin,side,matrix->first_row,\
			      matrix->findx,stream);
}

/************************************************************************************************/

template <class T>
void RGF<T>::FirstStage(int *b_min,int *b_max,int NBlock)
{
    int info;
    int IB;
    int NR,NM,NP;
    int NS1,NS2;
    T ONE   = f_one();
    T ZERO  = f_zero();
    int *ipiv = new int[b_size];

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);
        
    NR = b_max[0]-b_min[0];
	
    //extract E-H11
    extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,0,0);
    //eye11=M_dev
    init_eye_on_dev(M_dev,NR,0);
    //compute (E-H11-SigmaRl)^{-1}=M1_dev^{-1}=M_dev
    tgesv_dev(NR,NR,M1_dev,NR,ipiv,M_dev,NR,1,&info);
    
    for(IB=1;IB<NBlock-1;IB++){

        NR  = b_max[IB]-b_min[IB];
	NM  = b_max[IB-1]-b_min[IB-1];
	NS1 = b_min[IB];
	NS2 = b_min[IB-1];

	//extract -Hii-1=M2_dev
	extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,-1,stream_c);

	//compute -Hii-1*gi-1=M2_dev*M[add2]=M1_dev
	tgemm_dev(cublas_handle,'N','N',NR,NM,NM,ONE,M2_dev,NR,&M_dev[NS2*b_size],NM, \
		  ZERO,M1_dev,NR);
	//compute Hii-1*gi-1*Hi-1i=M1_dev*M2_dev'=M_dev[add1]
	tgemm_dev(cublas_handle,'N','C',NR,NR,NM,ONE,M1_dev,NR,M2_dev,NR,ZERO, \
		  &M_dev[NS1*b_size],NR);

	//extract E-Hii=M1_dev
	extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IB,stream_c);
	    
	//compute E-Hii-Hii-1*gi-1*Hi-1i=M1_dev-M_dev[add1]=M1_dev
	taxpy_dev(cublas_handle,NR*NR,-ONE,&M_dev[NS1*b_size],1,M1_dev,1);
       
	//eyeii=M_dev[add1]
	init_eye_on_dev(&M_dev[NS1*b_size],NR,0);
	    	    	    
	//compute (E-Hii-Hii-1*gi-1*Hi-1i)^{-1}=M_dev[add1]
	tgesv_dev(NR,NR,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,1,&info);
    }

    cudaStreamSynchronize(stream_c);
    cudaStreamSynchronize(stream_m);
    
    cublasSetStream((cublasHandle_t)cublas_handle,0);

    delete[] ipiv;
}

/************************************************************************************************/

template <class T>
void RGF<T>::SecondStage(T *GR,int *b_min,int *b_max,int NBlock)
{

    int info;
    int IB;
    int NR,NM,NP;
    int NS1,NS2;
    int N1,NN;
    T ONE      = f_one();
    T ZERO     = f_zero();
    int *ipiv    = new int[b_size];

    IB           = NBlock-1;
    NM           = b_max[IB-1]-b_min[IB-1];
    NR           = b_max[IB]-b_min[IB];
    NS1          = b_min[IB];
    NS2          = b_min[IB-1];

    //extract E-Hii=M3_dev
    extract_diag(edge_i_dev,index_j_dev,nnz_dev,M3_dev,b_min,b_max,IB,stream_c);

    //extract -Hii-1=M2_dev
    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,-1,stream_c);

    //compute -Hii-1*gi-1=M2_dev*M_dev[add2]=M1_dev
    tgemm_dev(cublas_handle,'N','N',NR,NM,NM,ONE,M2_dev,NR,&M_dev[NS2*b_size],NM, \
	      ZERO,M1_dev,NR);

    //compute E-Hii-Hii-1*gi-1*Hi-1i=M3_dev-M1_dev*M2_dev'=M3_dev
    tgemm_dev(cublas_handle,'N','C',NR,NR,NM,-ONE,M1_dev,NR,M2_dev,NR,ONE,M3_dev,NR);

    //eyeii=M_dev[add1]
    init_eye_on_dev(&M_dev[NS1*b_size],NR,0);

    //compute gi=(E-Hii-Hii-1*gi-1*Hi-1i)^{-1}=M3_dev^{-1}=M_dev[add1]
    tgesv_dev(NR,NR,M3_dev,NR,ipiv,&M_dev[NS1*b_size],NR,1,&info);

    //copy data to device
    cudaMemcpy(&GR[b_min[IB]*b_size],&M_dev[NS1*b_size],NR*NR*sizeof(T),\
	       cudaMemcpyDeviceToHost);

    for(IB=(NBlock-2);IB>=0;IB--){

        NR  = b_max[IB]-b_min[IB];
	NP  = b_max[IB+1]-b_min[IB+1];
	NS1 = b_min[IB];
	NS2 = b_min[IB+1];
	    
	//extract -Hii+1=M2_dev
	extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,1,stream_c);

	//compute -Hii+1*Gi+1=M2_dev*M_dev[add2]=M1_dev
	tgemm_dev(cublas_handle,'N','N',NR,NP,NP,ONE,M2_dev,NR,&M_dev[NS2*b_size],NP, \
		  ZERO,M1_dev,NR);
	//compute -gi*Hii+1*Gi+1=M_dev[add]*M1_dev=M3_dev=-Gii+1
	tgemm_dev(cublas_handle,'N','N',NR,NP,NR,ONE,&M_dev[NS1*b_size],NR,M1_dev,NR, \
		  ZERO,M3_dev,NR);
	//compute gi*Hii+1*Gi+1*Hi+1i=M3_dev*M2_dev'=M1_dev
	tgemm_dev(cublas_handle,'N','C',NR,NR,NP,ONE,M3_dev,NR,M2_dev,NR,ZERO,M1_dev,NR);
	//compute gi*Hii+1*Gi+1*Hi+1i*gi=M1_dev*M_dev[add1]=M2_dev
	tgemm_dev(cublas_handle,'N','N',NR,NR,NR,ONE,M1_dev,NR,&M_dev[NS1*b_size],NR, \
		  ZERO,M2_dev,NR);

	//Gi=gi+gi*Hii+1*Gi+1*Hi+1i*gi M_dev[add1]=M_dev[add1]+M2_dev
	taxpy_dev(cublas_handle,NR*NR,ONE,M2_dev,1,&M_dev[NS1*b_size],1);
	    
	//copy data to host
	cudaMemcpy(&GR[b_min[IB]*b_size],&M_dev[NS1*b_size],NR*NR*sizeof(T),\
		   cudaMemcpyDeviceToHost);
    }

    delete[] ipiv;
}

/************************************************************************************************/

template <class T>
void RGF<T>::FirstStageFactor(int *b_min,int *b_max,int NBlock)
{
    int info;
    int IB;
    int NR,NM;
    T ONE      = f_one();

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);

//tpotrf_dev('L', NR, M1_dev ,NR ,&info)
//rj dpotrf MF[0,0]
//rj dtrsm RLTN MF[0,0] MF[1,0]

    NR = b_max[0]-b_min[0];

    tpotrf_dev('L', NR, &MF_dev[mf_block_index(0, 0)], mf_block_lda(0, 0), &info);
    //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(0, 0), mf_block_lda(0, 0));

    //rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
    if (matrix->nt > 1)
    {
       ttrsm_dev('R', 'L', 'T', 'N', NR+matrix->nd, NR, ONE, &MF_dev[mf_block_index(0, 0)], mf_block_lda(0, 0), &MF_dev[mf_block_index(1, 0)], mf_block_lda(1, 0));
    }
    
    //rj IB = 1; IB < NBlock; IB++
    for (IB = 1; IB < matrix->nt-1; IB++)
    {
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NT M[IB,IB-1] M[IB,IB-1]
       // todo rj: 3-last parameter ZERO in PARDISO
       tgemm_dev(cublas_handle, 'N', 'T', NR, NR, NM, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));
       if (matrix->nd > 0)
       {
          tgemm_dev(cublas_handle, 'N', 'T', matrix->nd, NR, NM, -ONE, &MF_dev[mf_dense_block_index(IB-1)], mf_dense_block_lda(IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB));
       }
       //printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c: %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dpotrf MF[IB,IB]
       tpotrf_dev('L', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
       //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
       ttrsm_dev('R', 'L', 'T', 'N', NR+matrix->nd, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &MF_dev[mf_block_index(IB+1, IB)], mf_block_lda(IB+1, IB));
       //printf("RJ: ttrsm NR: %d, a: %d, lda: %d, b: %d, ldb: %d\n", NR, mf_block_index(IB-1, IB-1), mf_block_lda(IB-1, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1));
    }

    if (matrix->nt > 1)
    {
       IB = matrix->nt-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NT M[IB,IB-1] M[IB,IB-1]
       // todo rj: 3-last parameter ZERO in PARDISO
       tgemm_dev(cublas_handle, 'N', 'T', NR, NR, NM, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));
       if (matrix->nd > 0)
       {
          tgemm_dev(cublas_handle, 'N', 'T', matrix->nd, NR, NM, -ONE, &MF_dev[mf_dense_block_index(IB-1)], mf_dense_block_lda(IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB));
       }
       //printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c: %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dpotrf MF[IB,IB]
       tpotrf_dev('L', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
       //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
       ttrsm_dev('R', 'L', 'T', 'N', matrix->nd, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB));
       //printf("RJ: ttrsm NR: %d, a: %d, lda: %d, b: %d, ldb: %d\n", NR, mf_block_index(IB-1, IB-1), mf_block_lda(IB-1, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1));
    }

    if (matrix->nd > 0)
    {
       IB = NBlock-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NT M[IB,IB-1] M[IB,IB-1]
       //rj NM is the same for all blocks 0..nt-1
       for (int i = 0; i < NBlock-1; i++)
       {
          tgemm_dev(cublas_handle, 'N', 'T', NR, NR, NM, -ONE, &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));
       }
       //printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c: %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dpotrf MF[IB,IB]
       tpotrf_dev('L', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
       //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(IB, IB), mf_block_lda(IB, IB));
    }

    cudaStreamSynchronize(stream_c);
    cudaStreamSynchronize(stream_m);
    
    cublasSetStream((cublasHandle_t)cublas_handle,0);
}

/************************************************************************************************/

template <class T>
void RGF<T>::SecondStageSolve(int *b_min,int *b_max,int NBlock,int nrhs)
{
    int info;
    int IB;
    int NR,NM,NP;
    T ONE      = f_one();

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);

    //Forward pass
    NR = b_max[0]-b_min[0];

    //rj dtrsm LLNN MF[0,0] rhs[0]
    ttrsm_dev('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(0, 0)], mf_block_lda(0, 0), &rhs_dev[b_min[0]], matrix->size);

    //rj IB = 1; IB < NBlock; IB++
    for (IB = 1; IB < matrix->nt; IB++)
    {
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NN M[IB,IB-1] rhs[IB-1] rhs[IB]
       tgemm_dev(cublas_handle, 'N', 'N', NR, nrhs, NM, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &rhs_dev[b_min[IB-1]], matrix->size, ONE, &rhs_dev[b_min[IB]], matrix->size);
       //rj dtrsm LLNN MF[IB,IB] rhs[IB]
       ttrsm_dev('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size);
    }

    if (matrix->nd > 0)
    {
       IB = NBlock-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NN M[IB,IB-1] rhs[IB-1] rhs[IB]
       for (int i = 0; i < NBlock-1; i++)
       {
          tgemm_dev(cublas_handle, 'N', 'N', NR, nrhs, NM, -ONE, &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), &rhs_dev[b_min[i]], matrix->size, ONE, &rhs_dev[b_min[IB]], matrix->size);
       }
       //rj dtrsm LLNN MF[IB,IB] rhs[IB]
       ttrsm_dev('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size);
    }

    //Backward pass
    if (matrix->nd > 0)
    {
       IB = NBlock-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dtrsm LLTN MF[NBlock-1,NBlock-1] rhs[NBlock-1]
       ttrsm_dev('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size);

       //rj dgemm TN M[IB,IB-1] rhs[IB] rhs[IB-1]
       for (int i = 0; i < NBlock-1; i++)
       {
          tgemm_dev(cublas_handle, 'T', 'N', NM, nrhs, NR, -ONE, &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), &rhs_dev[b_min[IB]], matrix->size, ONE, &rhs_dev[b_min[i]], matrix->size);
       }
    }

    for (IB = matrix->nt-1; IB > 0; IB--)
    {
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dtrsm LLTN MF[IB,IB] rhs[IB]
       ttrsm_dev('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size);
       //rj dgemm TN M[IB,IB-1] rhs[IB] rhs[IB-1]
       tgemm_dev(cublas_handle, 'T', 'N', NM, nrhs, NR, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &rhs_dev[b_min[IB]], matrix->size, ONE, &rhs_dev[b_min[IB-1]], matrix->size);
    }

    IB = 0;
    NR = b_max[IB]-b_min[IB];

    //rj dtrsm LLTN MF[NBlock-1,NBlock-1] rhs[NBlock-1]
    ttrsm_dev('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size);

    cudaStreamSynchronize(stream_c);
    cudaStreamSynchronize(stream_m);
    
    cublasSetStream((cublasHandle_t)cublas_handle,0);
}

/************************************************************************************************/

template <class T>
void RGF<T>::ThirdStageRGF(T *tmp1_dev,T *tmp2_dev,int *b_min,int *b_max,int NBlock)
{
   int info;
   int IB;
   int NR,NM,NP;
   T ONE      = f_one();
   T ZERO     = f_zero();

   cublasSetStream((cublasHandle_t)cublas_handle,stream_c);

   if (matrix->nd > 0)
   {
      //dense block
      IB = NBlock-1;
      NR = b_max[IB]-b_min[IB];

      tril_dev(&MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), NR);
      ttrtri_dev('L', 'N', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
      tgemm_dev(cublas_handle, 'T', 'N', NR, NR, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), ZERO, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));

      //last non-dense block
      IB = NBlock-2;
      NR = b_max[IB]-b_min[IB];
      NP = b_max[IB+1]-b_min[IB+1];

      tril_dev(&MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), NR);
      ttrtri_dev('L', 'N', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);

      init_eye_on_dev(tmp1_dev,NR,0);

      tgemm_dev(cublas_handle, 'T', 'N', NR, NP, NP, ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), &MF_dev[mf_block_index(IB+1, IB+1)], mf_block_lda(IB+1, IB+1), ZERO, tmp2_dev, NR);
      tgemm_dev(cublas_handle, 'N', 'N', NR, NR, NP, ONE, tmp2_dev, NR, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), ONE, tmp1_dev, NR);

      tgemm_dev(cublas_handle, 'T', 'N', NR, NR, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), tmp1_dev, NR, ZERO, tmp2_dev, NR);
      tgemm_dev(cublas_handle, 'N', 'N', NR, NR, NR, ONE, tmp2_dev, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), ZERO, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));
   }
   else
   {
      //last non-dense block
      IB = NBlock-1;
      NR = b_max[IB]-b_min[IB];

      //memcpy_to_device(&aa,MF_dev+65,sizeof(double));
      tril_dev(&MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), NR);
      ttrtri_dev('L', 'N', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
      tgemm_dev(cublas_handle, 'T', 'N', NR, NR, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), ZERO, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));
   }

   //second-last non-dense block .. 0-block
   for (IB = matrix->nt-2; IB > -1; IB--)
   {
      NR = b_max[IB]-b_min[IB];
      NP = b_max[IB+1]-b_min[IB+1];

      tril_dev(&MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), NR);
      ttrtri_dev('L', 'N', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);

      init_eye_on_dev(tmp1_dev,NR,0);

      tgemm_dev(cublas_handle, 'T', 'N', NR, NP, NP, ONE, &MF_dev[mf_block_index(IB+1, IB)], mf_block_lda(IB+1, IB), &MF_dev[mf_block_index(IB+1, IB+1)], mf_block_lda(IB+1, IB+1), ZERO, tmp2_dev, NR);
      tgemm_dev(cublas_handle, 'N', 'N', NR, NR, NP, ONE, tmp2_dev, NR, &MF_dev[mf_block_index(IB+1, IB)], mf_block_lda(IB+1, IB), ONE, tmp1_dev, NR);

      tgemm_dev(cublas_handle, 'T', 'N', NR, NP, NP, ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), &MF_dev[mf_block_index(NBlock, NBlock)], mf_block_lda(NBlock, NBlock), ZERO, tmp2_dev, NR);
      tgemm_dev(cublas_handle, 'N', 'N', NR, NR, NP, ONE, tmp2_dev, NR, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), ONE, tmp1_dev, NR);

      tgemm_dev(cublas_handle, 'T', 'N', NR, NR, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), tmp1_dev, NR, ZERO, tmp2_dev, NR);
      tgemm_dev(cublas_handle, 'N', 'N', NR, NR, NR, ONE, tmp2_dev, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), ZERO, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB));
   }

   cudaStreamSynchronize(stream_c);
   cudaStreamSynchronize(stream_m);

   cublasSetStream((cublasHandle_t)cublas_handle,0);
}

/************************************************************************************************/

template <class T>
void RGF<T>::solve_equation(T *GR)
{
    int NBlock = matrix->NBlock;
    int *Bmin = matrix->Bmin;
    int *Bmax = matrix->Bmax;

    cudaStreamCreate(&stream_c);
    cudaStreamCreate(&stream_m);

    magma_init();
    cublas_init(&cublas_handle);

    create_blocks(Bmin,Bmax,NBlock);
    
    //Data allocation
    allocate_data_on_device((void**)&M_dev,b_size*matrix->size*sizeof(T));
    allocate_data_on_device((void**)&M1_dev,b_size*b_size*sizeof(T));
    allocate_data_on_device((void**)&M2_dev,b_size*b_size*sizeof(T));
    allocate_data_on_device((void**)&M3_dev,b_size*b_size*sizeof(T));

    allocate_data_on_device((void**)&edge_i_dev,(matrix->size+1)*sizeof(int));
    allocate_data_on_device((void**)&index_j_dev,matrix->n_nonzeros*sizeof(int));
    allocate_data_on_device((void**)&nnz_dev,matrix->n_nonzeros*sizeof(T));

    memcpy_to_device(matrix->edge_i,edge_i_dev,(matrix->size+1)*sizeof(int));
    memcpy_to_device(matrix->index_j,index_j_dev,matrix->n_nonzeros*sizeof(int));
    memcpy_to_device(matrix->nnz,nnz_dev,matrix->n_nonzeros*sizeof(T));

    //Computation
    init_var(GR,Bmax[NBlock-1]*b_size);
    
    FirstStage(Bmin,Bmax,NBlock);
    SecondStage(GR,Bmin,Bmax,NBlock);

    //Data deallocation
    deallocate_data_on_dev(edge_i_dev,(matrix->size+1)*sizeof(int));
    deallocate_data_on_dev(index_j_dev,matrix->n_nonzeros*sizeof(int));
    deallocate_data_on_dev(nnz_dev,matrix->n_nonzeros*sizeof(T));

    deallocate_data_on_dev(M_dev,b_size*matrix->size*sizeof(T));
    deallocate_data_on_dev(M1_dev,b_size*b_size*sizeof(T));
    deallocate_data_on_dev(M2_dev,b_size*b_size*sizeof(T));
    deallocate_data_on_dev(M3_dev,b_size*b_size*sizeof(T));

    magma_finalize();
    cublas_finalize(cublas_handle);

    cudaStreamDestroy(stream_c);
    cudaStreamDestroy(stream_m);
}

/************************************************************************************************/

template <class T>
void RGF<T>::factorize()
{
    int NBlock = matrix->NBlock;
    int *Bmin = matrix->Bmin;
    int *Bmax = matrix->Bmax;

    cudaStreamCreate(&stream_c);
    cudaStreamCreate(&stream_m);

    magma_init();
    cublas_init(&cublas_handle);
    //magma_queue_create_from_cuda(rj_device, stream_c, cublas_handle, NULL, &magma_queue);

    create_blocks(Bmin,Bmax,NBlock);
    
    //Data allocation
    allocate_data_on_device((void**)&MF_dev,matrix->n_nonzeros*sizeof(T));
    MF_dev_allocated = true;

    //Copy data to device
    memcpy_to_device(matrix->nnz,MF_dev,matrix->n_nonzeros*sizeof(T));

    //Computation
    FirstStageFactor(Bmin,Bmax,NBlock);

    //magma_queue_destroy(magma_queue);
    magma_finalize();
    cublas_finalize(cublas_handle);

    cudaStreamDestroy(stream_c);
    cudaStreamDestroy(stream_m);
}

/************************************************************************************************/

template <class T>
void RGF<T>::solve(T *b, int nrhs)
{
    solve(b, b, nrhs);
}

/************************************************************************************************/

template <class T>
void RGF<T>::solve(T *x, T *b, int nrhs)
{
    int NBlock = matrix->NBlock;
    int *Bmin = matrix->Bmin;
    int *Bmax = matrix->Bmax;

    cudaStreamCreate(&stream_c);
    cudaStreamCreate(&stream_m);

    magma_init();
    cublas_init(&cublas_handle);
    
    //Data allocation
    allocate_data_on_device((void**)&rhs_dev,matrix->size*nrhs*sizeof(T));

    //Copy data to device
    memcpy_to_device(b,rhs_dev,matrix->size*nrhs*sizeof(T));

    //Computation
    SecondStageSolve(Bmin,Bmax,NBlock,nrhs);

    //Copy data to host
    memcpy_to_host(x,rhs_dev,matrix->size*nrhs*sizeof(T));

    //Data deallocation
    deallocate_data_on_dev(rhs_dev,matrix->size*nrhs*sizeof(T));

    magma_finalize();
    cublas_finalize(cublas_handle);

    cudaStreamDestroy(stream_c);
    cudaStreamDestroy(stream_m);
}

/************************************************************************************************/

template <class T>
void RGF<T>::RGFdiag(T *diag)
{
    int NBlock = matrix->NBlock;
    int *Bmin = matrix->Bmin;
    int *Bmax = matrix->Bmax;

    cudaStreamCreate(&stream_c);
    cudaStreamCreate(&stream_m);

    magma_init();
    cublas_init(&cublas_handle);

    //Data allocation
    T *diag_dev;
    int *diag_pos_dev;
    T *tmp1_dev;
    T *tmp2_dev;
    allocate_data_on_device((void**)&diag_dev,matrix->size*sizeof(T));
    allocate_data_on_device((void**)&diag_pos_dev,matrix->size*sizeof(int));
    allocate_data_on_device((void**)&tmp1_dev,b_size*b_size*sizeof(T));
    allocate_data_on_device((void**)&tmp2_dev,b_size*b_size*sizeof(T));

    //Copy diag_pos to device
    memcpy_to_device(matrix->diag_pos,diag_pos_dev,matrix->size*sizeof(int));

    //Computation
    ThirdStageRGF(tmp1_dev,tmp2_dev,Bmin,Bmax,NBlock);
    indexed_copy_dev(MF_dev, diag_dev, diag_pos_dev, matrix->size);

    //Copy data to host
    memcpy_to_host(diag,diag_dev,matrix->size*sizeof(T));

    //Data deallocation
    deallocate_data_on_dev(diag_dev,matrix->size*sizeof(T));
    deallocate_data_on_dev(diag_pos_dev,matrix->size*sizeof(int));
    deallocate_data_on_dev(tmp1_dev,b_size*b_size*sizeof(T));
    deallocate_data_on_dev(tmp2_dev,b_size*b_size*sizeof(T));

    magma_finalize();
    cublas_finalize(cublas_handle);

    cudaStreamDestroy(stream_c);
    cudaStreamDestroy(stream_m);
}

/************************************************************************************************/

template <class T>
void RGF<T>::create_blocks(int *Bmin,int *Bmax,int NBlock)
{
    int IB;

    b_size = 0;

    for(IB=0;IB<NBlock;IB++){

        if(Bmax[IB]-Bmin[IB]>b_size){
	    b_size = Bmax[IB]-Bmin[IB];
	}
    }
}

/************************************************************************************************/

template <class T>
void RGF<T>::write_matrix(const char *filename,T *matrix,int NR,int NC)
{
    int IC,IR;
    ofstream myfile;
    
    myfile.open(filename);
    myfile.precision(8);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
            myfile<<real(matrix[IR+IC*NR])<<" "<<imag(matrix[IR+IC*NR])<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/

template <>
CPX RGF<CPX>::f_one()
{
    return CPX(1.0, 0.0);
}

template <>
double RGF<double>::f_one()
{
    return 1.0;
}

/************************************************************************************************/

template <>
CPX RGF<CPX>::f_zero()
{
    return CPX(0.0, 0.0);
}

template <>
double RGF<double>::f_zero()
{
    return 0.0;
}

/************************************************************************************************/

#endif
