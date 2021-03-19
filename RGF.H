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
    void factorize(size_t*,size_t*,T*);
    void solve(T*,size_t);
    void solve(T*,T*,size_t);
    void RGFdiag(T*);
    T logDet();
    double residualNorm(T*,T*);
    double residualNormNormalized(T*,T*);
				
private:

    TCSR<T>* matrix;
    int findx;
    size_t b_size;
    bool MF_dev_allocated = false;
    bool factorization_completed = false;

    magma_queue_t magma_queue;
    cudaStream_t stream_c;
    void *cublas_handle;

    T *M_dev;
    T *M1_dev;
    T *M2_dev;
    T *M3_dev;
    T *MF_dev;
    T *rhs_dev;
    T *nnz_dev;
    size_t *edge_i_dev;
    size_t *index_j_dev;

    inline size_t mf_block_index(size_t, size_t);
    inline size_t mf_block_lda(size_t, size_t);
    inline size_t mf_dense_block_index(size_t);
    inline size_t mf_dense_block_lda(size_t);

    void extract_diag(int*,int*,T*,T*,int*,int*,int,cudaStream_t);
    void extract_not_diag(int*,int*,T*,T*,int*,int*,int,int,cudaStream_t);
    void FirstStage(int*,int*,int);
    void SecondStage(T*,int*,int*,int);
    void FirstStageFactor(size_t*,size_t*,size_t);
    void SecondStageSolve(size_t*,size_t*,size_t,size_t);
    void ThirdStageRGF(T*,T*,size_t*,size_t*,size_t);
    void create_blocks(size_t*,size_t*,size_t);
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

    magma_init();
    magma_device_t device;
    magma_getdevice(&device);
    magma_queue_create_from_cuda(device, NULL, NULL, NULL, &magma_queue);
    stream_c = magma_queue_get_cuda_stream(magma_queue);
    cublas_handle = magma_queue_get_cublas_handle(magma_queue);
}

/************************************************************************************************/

template <class T>
RGF<T>::~RGF()
{
   if (MF_dev_allocated)
   {
      deallocate_data_on_dev(MF_dev,matrix->n_nonzeros*sizeof(T));
   }

   magma_queue_destroy(magma_queue);
   magma_finalize();
}

/************************************************************************************************/

template <class T>
inline size_t RGF<T>::mf_block_index(size_t r, size_t c)
{
   return matrix->diag_pos[c*b_size] + (r-c)*b_size;
}

/************************************************************************************************/

template <class T>
inline size_t RGF<T>::mf_block_lda(size_t r, size_t c)
{
   return matrix->index_i[c*b_size];
}

/************************************************************************************************/

template <class T>
inline size_t RGF<T>::mf_dense_block_index(size_t i)
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
inline size_t RGF<T>::mf_dense_block_lda(size_t i)
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
	tgemm_dev('N','N',NR,NM,NM,ONE,M2_dev,NR,&M_dev[NS2*b_size],NM, \
		  ZERO,M1_dev,NR,magma_queue);
	//compute Hii-1*gi-1*Hi-1i=M1_dev*M2_dev'=M_dev[add1]
	tgemm_dev('N','C',NR,NR,NM,ONE,M1_dev,NR,M2_dev,NR,ZERO, \
		  &M_dev[NS1*b_size],NR,magma_queue);

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
    magma_queue_sync(magma_queue);
    
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
    tgemm_dev('N','N',NR,NM,NM,ONE,M2_dev,NR,&M_dev[NS2*b_size],NM, \
	      ZERO,M1_dev,NR,magma_queue);

    //compute E-Hii-Hii-1*gi-1*Hi-1i=M3_dev-M1_dev*M2_dev'=M3_dev
    tgemm_dev('N','C',NR,NR,NM,-ONE,M1_dev,NR,M2_dev,NR,ONE,M3_dev,NR,magma_queue);

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
	tgemm_dev('N','N',NR,NP,NP,ONE,M2_dev,NR,&M_dev[NS2*b_size],NP, \
		  ZERO,M1_dev,NR,magma_queue);
	//compute -gi*Hii+1*Gi+1=M_dev[add]*M1_dev=M3_dev=-Gii+1
	tgemm_dev('N','N',NR,NP,NR,ONE,&M_dev[NS1*b_size],NR,M1_dev,NR, \
		  ZERO,M3_dev,NR,magma_queue);
	//compute gi*Hii+1*Gi+1*Hi+1i=M3_dev*M2_dev'=M1_dev
	tgemm_dev('N','C',NR,NR,NP,ONE,M3_dev,NR,M2_dev,NR,ZERO,M1_dev,NR,magma_queue);
	//compute gi*Hii+1*Gi+1*Hi+1i*gi=M1_dev*M_dev[add1]=M2_dev
	tgemm_dev('N','N',NR,NR,NR,ONE,M1_dev,NR,&M_dev[NS1*b_size],NR, \
		  ZERO,M2_dev,NR,magma_queue);

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
void RGF<T>::FirstStageFactor(size_t *b_min,size_t *b_max,size_t NBlock)
{
    int info;
    size_t IB;
    size_t NR,NM;
    T ONE      = f_one();

//tpotrf_dev('L', NR, M1_dev ,NR ,&info)
//rj dpotrf MF[0,0]
//rj dtrsm RLTN MF[0,0] MF[1,0]

    NR = b_max[0]-b_min[0];

    tpotrf_dev('L', NR, &MF_dev[mf_block_index(0, 0)], mf_block_lda(0, 0), &info);
    //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(0, 0), mf_block_lda(0, 0));

    //rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
    if (matrix->nt > 1)
    {
       ttrsm_dev('R', 'L', 'T', 'N', NR+matrix->nd, NR, ONE, &MF_dev[mf_block_index(0, 0)], mf_block_lda(0, 0), &MF_dev[mf_block_index(1, 0)], mf_block_lda(1, 0), magma_queue);
    }
    
    //rj IB = 1; IB < NBlock; IB++
    for (IB = 1; IB < matrix->nt-1; IB++)
    {
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NT M[IB,IB-1] M[IB,IB-1]
       // todo rj: 3-last parameter ZERO in PARDISO
       tgemm_dev('N', 'T', NR, NR, NM, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB),magma_queue);
       if (matrix->nd > 0)
       {
          tgemm_dev('N', 'T', matrix->nd, NR, NM, -ONE, &MF_dev[mf_dense_block_index(IB-1)], mf_dense_block_lda(IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB),magma_queue);
       }
       //printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c: %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dpotrf MF[IB,IB]
       tpotrf_dev('L', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
       //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
       ttrsm_dev('R', 'L', 'T', 'N', NR+matrix->nd, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &MF_dev[mf_block_index(IB+1, IB)], mf_block_lda(IB+1, IB), magma_queue);
       //printf("RJ: ttrsm NR: %d, a: %d, lda: %d, b: %d, ldb: %d\n", NR, mf_block_index(IB-1, IB-1), mf_block_lda(IB-1, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1));
    }

    if (matrix->nt > 1)
    {
       IB = matrix->nt-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NT M[IB,IB-1] M[IB,IB-1]
       // todo rj: 3-last parameter ZERO in PARDISO
       tgemm_dev('N', 'T', NR, NR, NM, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB),magma_queue);
       if (matrix->nd > 0)
       {
          tgemm_dev('N', 'T', matrix->nd, NR, NM, -ONE, &MF_dev[mf_dense_block_index(IB-1)], mf_dense_block_lda(IB-1), &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB),magma_queue);
       }
       //printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c: %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dpotrf MF[IB,IB]
       tpotrf_dev('L', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
       //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
       ttrsm_dev('R', 'L', 'T', 'N', matrix->nd, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), magma_queue);
       //printf("RJ: ttrsm NR: %d, a: %d, lda: %d, b: %d, ldb: %d\n", NR, mf_block_index(IB-1, IB-1), mf_block_lda(IB-1, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1));
    }

    if (matrix->nd > 0)
    {
       IB = NBlock-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NT M[IB,IB-1] M[IB,IB-1]
       //rj NM is the same for all blocks 0..nt-1
       for (size_t i = 0; i < NBlock-1; i++)
       {
          tgemm_dev('N', 'T', NR, NR, NM, -ONE, &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB),magma_queue);
       }
       //printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c: %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1), mf_block_index(IB, IB), mf_block_lda(IB, IB));
       //rj dpotrf MF[IB,IB]
       tpotrf_dev('L', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);
       //printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(IB, IB), mf_block_lda(IB, IB));
    }

    cudaStreamSynchronize(stream_c);
    magma_queue_sync(magma_queue);
}

/************************************************************************************************/

template <class T>
void RGF<T>::SecondStageSolve(size_t *b_min,size_t *b_max,size_t NBlock,size_t nrhs)
{
    int info;
    size_t IB;
    size_t NR,NM,NP;
    T ONE      = f_one();

    //Forward pass
    NR = b_max[0]-b_min[0];

    //rj dtrsm LLNN MF[0,0] rhs[0]
    ttrsm_dev('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(0, 0)], mf_block_lda(0, 0), &rhs_dev[b_min[0]], matrix->size, magma_queue);

    //rj IB = 1; IB < NBlock; IB++
    for (IB = 1; IB < matrix->nt; IB++)
    {
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NN M[IB,IB-1] rhs[IB-1] rhs[IB]
       tgemm_dev('N', 'N', NR, nrhs, NM, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &rhs_dev[b_min[IB-1]], matrix->size, ONE, &rhs_dev[b_min[IB]], matrix->size,magma_queue);
       //rj dtrsm LLNN MF[IB,IB] rhs[IB]
       ttrsm_dev('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size, magma_queue);
    }

    if (matrix->nd > 0)
    {
       IB = NBlock-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dgemm NN M[IB,IB-1] rhs[IB-1] rhs[IB]
       for (size_t i = 0; i < NBlock-1; i++)
       {
          tgemm_dev('N', 'N', NR, nrhs, NM, -ONE, &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), &rhs_dev[b_min[i]], matrix->size, ONE, &rhs_dev[b_min[IB]], matrix->size,magma_queue);
       }
       //rj dtrsm LLNN MF[IB,IB] rhs[IB]
       ttrsm_dev('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size, magma_queue);
    }

    //Backward pass
    if (matrix->nd > 0)
    {
       IB = NBlock-1;
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dtrsm LLTN MF[NBlock-1,NBlock-1] rhs[NBlock-1]
       ttrsm_dev('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size, magma_queue);

       //rj dgemm TN M[IB,IB-1] rhs[IB] rhs[IB-1]
       for (size_t i = 0; i < NBlock-1; i++)
       {
          tgemm_dev('T', 'N', NM, nrhs, NR, -ONE, &MF_dev[mf_dense_block_index(i)], mf_dense_block_lda(i), &rhs_dev[b_min[IB]], matrix->size, ONE, &rhs_dev[b_min[i]], matrix->size,magma_queue);
       }
    }

    for (IB = matrix->nt-1; IB > 0; IB--)
    {
       NR  = b_max[IB]-b_min[IB];
       NM  = b_max[IB-1]-b_min[IB-1];

       //rj dtrsm LLTN MF[IB,IB] rhs[IB]
       ttrsm_dev('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size, magma_queue);
       //rj dgemm TN M[IB,IB-1] rhs[IB] rhs[IB-1]
       tgemm_dev('T', 'N', NM, nrhs, NR, -ONE, &MF_dev[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &rhs_dev[b_min[IB]], matrix->size, ONE, &rhs_dev[b_min[IB-1]], matrix->size,magma_queue);
    }

    IB = 0;
    NR = b_max[IB]-b_min[IB];

    //rj dtrsm LLTN MF[NBlock-1,NBlock-1] rhs[NBlock-1]
    ttrsm_dev('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &rhs_dev[b_min[IB]], matrix->size, magma_queue);

    cudaStreamSynchronize(stream_c);
    magma_queue_sync(magma_queue);
}

/************************************************************************************************/

template <class T>
void RGF<T>::ThirdStageRGF(T *tmp1_dev,T *tmp2_dev,size_t *b_min,size_t *b_max,size_t NBlock)
{
   int info;
   size_t IB;
   size_t NR,NM,NP,ND;
   T ONE      = f_one();
   T ZERO     = f_zero();

   ND = matrix->nd;

   if (matrix->nd > 0)
   {
      //dense block
      IB = NBlock-1;
      NR = b_max[IB]-b_min[IB];

      tlacpy_dev('L', NR, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), tmp1_dev, NR, magma_queue);
      tril_dev(tmp1_dev, NR, NR);
      ttrtri_dev('L', 'N', NR, tmp1_dev, NR, &info);
      tgemm_dev('T', 'N', NR, NR, NR, ONE, tmp1_dev, NR, tmp1_dev, NR, ZERO, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB),magma_queue);

      //last non-dense block
      IB = NBlock-2;
      NR = b_max[IB]-b_min[IB];
      NP = b_max[IB+1]-b_min[IB+1];

      tril_dev(&MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), NR);
      ttrtri_dev('L', 'N', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);

      init_eye_on_dev(tmp1_dev,NR,0);

      tgemm_dev('T', 'N', NR, NP, NP, ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), &MF_dev[mf_block_index(IB+1, IB+1)], mf_block_lda(IB+1, IB+1), ZERO, tmp2_dev, NR,magma_queue);
      tgemm_dev('N', 'N', NR, NR, NP, ONE, tmp2_dev, NR, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), ONE, tmp1_dev, NR,magma_queue);

      tgemm_dev('T', 'N', NR, NR, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), tmp1_dev, NR, ZERO, tmp2_dev, NR,magma_queue);
      tgemm_dev('N', 'N', NR, NR, NR, ONE, tmp2_dev, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), ZERO, tmp1_dev, NR,magma_queue);
      tlacpy_dev('F', NR, NR, tmp1_dev, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), magma_queue);
   }
   else
   {
      //last non-dense block
      IB = NBlock-1;
      NR = b_max[IB]-b_min[IB];

      tlacpy_dev('L', NR, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), tmp1_dev, NR, magma_queue);
      tril_dev(tmp1_dev, NR, NR);
      ttrtri_dev('L', 'N', NR, tmp1_dev, NR, &info);
      tgemm_dev('T', 'N', NR, NR, NR, ONE, tmp1_dev, NR, tmp1_dev, NR, ZERO, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB),magma_queue);
   }

   //second-last non-dense block .. 0-block
   for (int IBi = matrix->nt-2; IBi > -1; IBi--)
   {
      IB = IBi;
      NR = b_max[IB]-b_min[IB];
      NP = b_max[IB+1]-b_min[IB+1];

      tril_dev(&MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), NR);
      ttrtri_dev('L', 'N', NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &info);

      init_eye_on_dev(tmp1_dev,NR,0);

      tgemm_dev('T', 'N', NR, NP, NP, ONE, &MF_dev[mf_block_index(IB+1, IB)], mf_block_lda(IB+1, IB), &MF_dev[mf_block_index(IB+1, IB+1)], mf_block_lda(IB+1, IB+1), ZERO, tmp2_dev, NR,magma_queue);
      tgemm_dev('N', 'N', NR, NR, NP, ONE, tmp2_dev, NR, &MF_dev[mf_block_index(IB+1, IB)], mf_block_lda(IB+1, IB), ONE, tmp1_dev, NR,magma_queue);

      if (matrix->nd > 0)
      {
         tgemm_dev('T', 'N', NR, ND, ND, ONE, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), &MF_dev[mf_block_index(NBlock-1, NBlock-1)], mf_block_lda(NBlock-1, NBlock-1), ZERO, tmp2_dev, NR,magma_queue);
         tgemm_dev('N', 'N', NR, NR, ND, ONE, tmp2_dev, NR, &MF_dev[mf_dense_block_index(IB)], mf_dense_block_lda(IB), ONE, tmp1_dev, NR,magma_queue);
      }

      tgemm_dev('T', 'N', NR, NR, NR, ONE, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), tmp1_dev, NR, ZERO, tmp2_dev, NR,magma_queue);
      tgemm_dev('N', 'N', NR, NR, NR, ONE, tmp2_dev, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), ZERO, tmp1_dev, NR,magma_queue);
      tlacpy_dev('F', NR, NR, tmp1_dev, NR, &MF_dev[mf_block_index(IB, IB)], mf_block_lda(IB, IB), magma_queue);
   }

   cudaStreamSynchronize(stream_c);
    magma_queue_sync(magma_queue);
}

/************************************************************************************************/

template <class T>
void RGF<T>::solve_equation(T *GR)
{
    int NBlock = matrix->NBlock;
    int *Bmin = matrix->Bmin;
    int *Bmax = matrix->Bmax;

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
}

/************************************************************************************************/

template <class T>
void RGF<T>::factorize()
{
    size_t NBlock = matrix->NBlock;
    size_t *Bmin = matrix->Bmin;
    size_t *Bmax = matrix->Bmax;

    create_blocks(Bmin,Bmax,NBlock);

    //Data allocation
    allocate_data_on_device((void**)&MF_dev,matrix->n_nonzeros*sizeof(T));
    MF_dev_allocated = true;

    //Copy data to device
    memcpy_to_device(matrix->nnz,MF_dev,matrix->n_nonzeros*sizeof(T));

    //Computation
    FirstStageFactor(Bmin,Bmax,NBlock);

    factorization_completed = true;
}

/************************************************************************************************/

template <class T>
void RGF<T>::factorize(size_t *ia, size_t *ja, T *a)
{
    size_t NBlock = matrix->NBlock;
    size_t *Bmin = matrix->Bmin;
    size_t *Bmax = matrix->Bmax;

    create_blocks(Bmin,Bmax,NBlock);

    //Data allocation
    allocate_data_on_device((void**)&MF_dev,matrix->n_nonzeros*sizeof(T));
    MF_dev_allocated = true;
    size_t nnz = ia[matrix->size];
    size_t *ia_dev;
    size_t *ja_dev;
    T *a_dev;
    allocate_data_on_device((void**)&ia_dev,(matrix->size+1)*sizeof(size_t));
    allocate_data_on_device((void**)&ja_dev,nnz*sizeof(size_t));
    allocate_data_on_device((void**)&a_dev,nnz*sizeof(T));

    //Copy data to device
    memcpy_to_device(ia,ia_dev,(matrix->size+1)*sizeof(size_t));
    memcpy_to_device(ja,ja_dev,nnz*sizeof(size_t));
    memcpy_to_device(a,a_dev,nnz*sizeof(T));

    //Initialize block matrix on device
    fill_dev(MF_dev, T(0.0), matrix->n_nonzeros);
    init_block_matrix_dev(MF_dev, ia_dev, ja_dev, a_dev, nnz, matrix->ns, matrix->nt, matrix->nd);

    //Data deallocation
    deallocate_data_on_dev(ia_dev,(matrix->size+1)*sizeof(size_t));
    deallocate_data_on_dev(ja_dev,nnz*sizeof(size_t));
    deallocate_data_on_dev(a_dev,nnz*sizeof(T));

    //Computation
    FirstStageFactor(Bmin,Bmax,NBlock);

    factorization_completed = true;
}

/************************************************************************************************/

template <class T>
void RGF<T>::solve(T *b, size_t nrhs)
{
    solve(b, b, nrhs);
}

/************************************************************************************************/

template <class T>
void RGF<T>::solve(T *x, T *b, size_t nrhs)
{
    if (!factorization_completed)
       factorize();

    size_t NBlock = matrix->NBlock;
    size_t *Bmin = matrix->Bmin;
    size_t *Bmax = matrix->Bmax;
    
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
}

/************************************************************************************************/

template <class T>
void RGF<T>::RGFdiag(T *diag)
{
    if (!factorization_completed)
       factorize();

    factorization_completed = false;

    size_t NBlock = matrix->NBlock;
    size_t *Bmin = matrix->Bmin;
    size_t *Bmax = matrix->Bmax;

    //Data allocation
    T *diag_dev;
    size_t *diag_pos_dev;
    T *tmp1_dev;
    T *tmp2_dev;
    allocate_data_on_device((void**)&diag_dev,matrix->size*sizeof(T));
    allocate_data_on_device((void**)&diag_pos_dev,matrix->size*sizeof(size_t));
    allocate_data_on_device((void**)&tmp1_dev,b_size*b_size*sizeof(T));
    allocate_data_on_device((void**)&tmp2_dev,b_size*b_size*sizeof(T));

    //Copy diag_pos to device
    memcpy_to_device(matrix->diag_pos,diag_pos_dev,matrix->size*sizeof(size_t));

    //Computation
    ThirdStageRGF(tmp1_dev,tmp2_dev,Bmin,Bmax,NBlock);
    indexed_copy_dev(MF_dev, diag_dev, diag_pos_dev, matrix->size);

    //Copy data to host
    memcpy_to_host(diag,diag_dev,matrix->size*sizeof(T));

    //Data deallocation
    deallocate_data_on_dev(diag_dev,matrix->size*sizeof(T));
    deallocate_data_on_dev(diag_pos_dev,matrix->size*sizeof(size_t));
    deallocate_data_on_dev(tmp1_dev,b_size*b_size*sizeof(T));
    deallocate_data_on_dev(tmp2_dev,b_size*b_size*sizeof(T));
}

/************************************************************************************************/

template <class T>
T RGF<T>::logDet()
{
    if (!factorization_completed)
       factorize();

    //Data allocation
    T *diag_dev;
    size_t *diag_pos_dev;
    allocate_data_on_device((void**)&diag_dev,matrix->size*sizeof(T));
    allocate_data_on_device((void**)&diag_pos_dev,matrix->size*sizeof(size_t));
    T det = T(0.0);

    //Copy diag_pos to device
    memcpy_to_device(matrix->diag_pos,diag_pos_dev,matrix->size*sizeof(size_t));

    //Computation
    indexed_copy_dev(MF_dev, diag_dev, diag_pos_dev, matrix->size);
    log_dev(diag_dev, matrix->size);
    tsum_dev(cublas_handle, matrix->size, diag_dev, 1, &det);

    //Copy data to host
    //memcpy_to_host(diag,diag_dev,matrix->size*sizeof(T));

    //Data deallocation
    deallocate_data_on_dev(diag_dev,matrix->size*sizeof(T));
    deallocate_data_on_dev(diag_pos_dev,matrix->size*sizeof(size_t));

    return 2.0*det;
}

/************************************************************************************************/

template <class T>
double RGF<T>::residualNorm(T *x, T *b)
{
   T *r = new T[matrix->size];
   size_t NR, NM;
   size_t IB;
   size_t *b_min = matrix->Bmin;
   size_t *b_max = matrix->Bmax;
   size_t NBlock = matrix->NBlock;

   memcpy(r, b, matrix->size*sizeof(T));

   IB = 0;
   NR = b_max[IB]-b_min[IB];

   c_tsymm('L', 'L', NR, 1, -1.0, &matrix->nnz[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &x[b_min[IB]], matrix->size, 1.0, &r[b_min[IB]], matrix->size);

   for (size_t IB = 1; IB < matrix->nt; IB++)
   {
      NR = b_max[IB]-b_min[IB];
      NM = b_max[IB-1]-b_min[IB-1];

      c_tgemm('N', 'N', NR, 1, NM, -1.0, &matrix->nnz[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &x[b_min[IB-1]], matrix->size, 1.0, &r[b_min[IB]], matrix->size);
      c_tgemm('T', 'N', NM, 1, NR, -1.0, &matrix->nnz[mf_block_index(IB, IB-1)], mf_block_lda(IB, IB-1), &x[b_min[IB]], matrix->size, 1.0, &r[b_min[IB-1]], matrix->size);
      c_tsymm('L', 'L', NR, 1, -1.0, &matrix->nnz[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &x[b_min[IB]], matrix->size, 1.0, &r[b_min[IB]], matrix->size);
   }

   if (matrix->nd > 0)
   {
      NR = matrix->nd;

      for (size_t IB = 0; IB < matrix->nt; IB++)
      {
         NM = b_max[IB]-b_min[IB];

         c_tgemm('N', 'N', NR, 1, NM, -1.0, &matrix->nnz[mf_dense_block_index(IB)], mf_dense_block_lda(IB), &x[b_min[IB]], matrix->size, 1.0, &r[b_min[NBlock-1]], matrix->size);
         c_tgemm('T', 'N', NM, 1, NR, -1.0, &matrix->nnz[mf_dense_block_index(IB)], mf_dense_block_lda(IB), &x[b_min[NBlock-1]], matrix->size, 1.0, &r[b_min[IB]], matrix->size);
      }

      IB = NBlock-1;
      c_tsymm('L', 'L', NR, 1, -1.0, &matrix->nnz[mf_block_index(IB, IB)], mf_block_lda(IB, IB), &x[b_min[IB]], matrix->size, 1.0, &r[b_min[IB]], matrix->size);
   }

   double res = c_dtnrm2(matrix->size, r, 1);

   delete[] r;

   return res;
}

/************************************************************************************************/

template <class T>
double RGF<T>::residualNormNormalized(T *x, T *b)
{
   return residualNorm(x, b) / c_dtnrm2(matrix->size, b, 1);
}

/************************************************************************************************/

template <class T>
void RGF<T>::create_blocks(size_t *Bmin,size_t *Bmax,size_t NBlock)
{
    size_t IB;

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
