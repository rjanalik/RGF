#ifndef __RGF
#define __RGF

#include <string.h>

#include "CWC_utility.H"
#include "Types.H"
#include "Utilities.H"
#include "cublas_v2.h"

template <class T> class RGF {
  /**
   * @brief Solves blockdiagonal (arrow) matrices using the RGF algorithm
   * @details Performs selective inversion using the RGF algorithm
   * \note: atm returns only inverse diagonal elements
   * @param[inout] diag Description
   */

public:
  RGF(size_t *, size_t *, T *, size_t, size_t, size_t);

  ~RGF();
  double factorize();
  double solve(T *, size_t);
  double solve(T *, T *, size_t);
  double RGFdiag(T *);
  T logDet();
  double residualNorm(T *, T *);
  double residualNormNormalized(T *, T *);
  double flop_count_factorise();

private:
  size_t *matrix_ia;
  size_t *matrix_ja;
  T *matrix_a;
  size_t matrix_size; /**< matrix_size: Total matrix size of size \p ns * \p
                         nt + \p nd*/
  size_t matrix_nnz;  /**< total number of non-zeros of the sparse matrix A*/
  size_t matrix_ns;   /**< same as ns*/
  size_t matrix_nt;   /**< same as nt*/
  size_t matrix_nd;   /**< same as nd*/
  size_t *Bmin;  /**< Bmin: Array of size NBlock with the first row start idx of
                    the blocks*/
  size_t *Bmax;  /**< Bmax: Array of size NBlock with the first row end idx of
                    the block*/
  size_t NBlock; /**< NBlock: Number of blocks in x/y with nt or nt+1 (with
                    random effects)*/
  size_t *diag_pos;

  size_t b_size;                 /**< b_size: Maximum possible blocksize*/
  bool MF_dev_allocated = false; /**< MF_dev_allocated: Flag if matrix
                                    factorization device is allocated*/
  bool factorization_completed = false; /**< factorization_completed: Flag if
                                           factorization is completed*/

  magma_queue_t magma_queue;
  cudaStream_t stream_c;
  void *cublas_handle;

  /**
   * @name Different blocks
   * Description: pointers of the blocks get swaped each round.
   * @{
   */
  T *MF;         /**< MF: matrix Factorizatoin Storage buffer */
  T *blockR_dev; /**< blockR_dev: first calculates A=LL^T and stores L^T; then
                    calculates and stores trms(L_i^T, B_i)*/
  T *blockR_dev1; /**< blockR_dev: first calculates A=LL^T and stores L^T; then
                    calculates and stores trms(L_i^T, B_i)*/
  T *blockR_dev2; /**< blockR_dev: first calculates A=LL^T and stores L^T; then
                    calculates and stores trms(L_i^T, B_i)*/
  T *blockM_dev; /**< blockM_dev: Maximum possible blocksize*/
  T *blockP_dev; /**< blockP_dev: Maximum possible blocksize*/
  T *blockDense_dev;
  T *rhs;
  /**
   * @}
   */
  size_t *ia_dev;
  size_t *ja_dev;
  T *a_dev;
  T *diag_dev;
  size_t *diag_pos_dev;

  inline size_t mf_block_index(size_t, size_t);
  inline size_t mf_block_lda(size_t, size_t);
  inline size_t mf_dense_block_index(size_t);
  inline size_t mf_dense_block_offset(size_t);
  inline size_t mf_dense_block_lda(size_t);

  double FirstStageFactor();
  double SecondStageSolve(size_t);
  double ThirdStageRGF(T *, T *);
  void create_blocks();
  inline void init_supernode(T *M_dev, size_t supernode);
  inline void init_supernode(T *M_dev, size_t supernode, cudaStream_t stream);
  inline void copy_supernode_to_host(T *M_dev, size_t supernode,
                                     cudaStream_t stream = NULL);
  inline void copy_supernode_to_device(T *M_dev, size_t supernode,
                                       cudaStream_t stream = NULL);
  inline void copy_supernode_diag(T *src, size_t supernode,
                                  cudaStream_t stream = NULL);
  inline void swap_pointers(T **ptr1, T **ptr2);

  T f_one();
  T f_zero();
};

/************************************************************************************************/

/**
 * @brief Recursive Green Function Algorithm
 * @details Description
 * \todo should CSR not be entered as (val, col_idx, row_ptr)?
 * CSR format of sparse matrix A
 * @param[inout] ia row_ptr of size(row_ptr)=n+1
 * @param[out] ja row_idx of size(row_ptr)=nnz(A)
 * @param[out] a values of size(a)=nnz(A)
 * @param[in] ns number of spatial grid points
 * @param[in] nt temporal grid points
 * @param[in] nd/nb number of columns/rows/fixed effects
 * @return instance of RGF class
 */
// size_t size = ns * nt + nd;
template <class T>
RGF<T>::RGF(size_t *ia, size_t *ja, T *a, size_t ns, size_t nt, size_t nd) {
#ifdef DEBUG
  print_header("RGF<T>::RGF(size_t *ia, size_t *ja, T *a, size_t ns, size_t "
               "nt, size_t nd)");
#endif
  matrix_ia = ia;
  matrix_ja = ja;
  matrix_a = a;
  matrix_ns = ns;
  matrix_nt = nt;
  matrix_nd = nd;

  matrix_size = ns * nt + nd;
  matrix_nnz = 2 * (nt - 1) * ns * ns + ns * ns + matrix_size * nd;

  // NOTE: total number of blocks
  if (nd > 0)
    NBlock = nt + 1;
  else
    NBlock = nt;
  Bmin = new size_t[NBlock];
  Bmax = new size_t[NBlock];
  // Calculate the start- and end index of first row of each block
  for (size_t i = 0; i < nt; i++) {
    Bmin[i] = i * ns;
    Bmax[i] = (i + 1) * ns;
  }
  if (nd > 0) {
    Bmin[nt] = nt * ns;
    Bmax[nt] = nt * ns + nd;
  }

  // TODO: what do we fill the blocks here with and what is diag_pos?
  // blocks with lda = 2*ns + nd
  diag_pos = new size_t[matrix_size];
  size_t IB; /** IB = Index Block */
  // Loop through all the blocks
  for (IB = 0; IB < nt - 1; IB++) {
    // For a given block step through all elements
    for (size_t i = 0; i < ns; i++) {
      diag_pos[IB * ns + i] = IB * ns * (2 * ns + nd) + i * (2 * ns + nd + 1);
    }
  }
  // last block with lda = ns + nd
  IB = nt - 1;
  for (size_t i = 0; i < ns; i++) {
    diag_pos[IB * ns + i] = IB * ns * (2 * ns + nd) + i * (ns + nd + 1);
  }
  // dense block with lda = ns + nd
  if (nd > 0) {
    IB = nt;
    for (size_t i = 0; i < nd; i++) {
      diag_pos[nt * ns + i] =
          (nt - 1) * ns * (2 * ns + nd) + ns * (ns + nd) + i * (nd + 1);
    }
  }

  magma_init();
  magma_device_t device;
  magma_getdevice(&device);
  // TODO: maybe create different queues already here?
  magma_queue_create_from_cuda(device, NULL, NULL, NULL, &magma_queue);
  stream_c = magma_queue_get_cuda_stream(magma_queue);
  cublas_handle = magma_queue_get_cublas_handle(magma_queue);
}

/************************************************************************************************/

template <class T> RGF<T>::~RGF() {
#ifdef DEBUG
  print_header("~RGF()");
#endif
  if (MF_dev_allocated) {
    /**< max_supernode_nnz: size of supernode blocks
     * for nt=0 we have just one single block
     * for nt>0 we have a big supernode block
     * TODO: is a supernode block [B A fixed]
     */
    size_t max_supernode_nnz = matrix_nt > 1
                                   ? matrix_ns * (2 * matrix_ns + matrix_nd)
                                   : matrix_ns * (matrix_ns + matrix_nd);
    /**< dense_supernode_nnz: size of dense random effects blocks*/
    size_t dense_supernode_nnz = matrix_nd > 0 ? matrix_nd * matrix_nd : 0;
    deallocate_data_on_dev(blockR_dev, max_supernode_nnz * sizeof(T));
    deallocate_data_on_dev(blockR_dev1, max_supernode_nnz * sizeof(T));
    deallocate_data_on_dev(blockR_dev2, max_supernode_nnz * sizeof(T));
    deallocate_data_on_dev(blockM_dev, max_supernode_nnz * sizeof(T));
    deallocate_data_on_dev(blockP_dev, max_supernode_nnz * sizeof(T));
    deallocate_data_on_dev(blockDense_dev, dense_supernode_nnz * sizeof(T));
    delete[] MF;
  }

  magma_queue_destroy(magma_queue);
  magma_finalize();

  delete[] Bmin;
  delete[] Bmax;
  delete[] diag_pos;
}

/************************************************************************************************/

template <class T> inline size_t RGF<T>::mf_block_index(size_t r, size_t c) {
  return diag_pos[c * b_size] + (r - c) * b_size;
}

/************************************************************************************************/

template <class T>
/**
 * @brief returns the correct leading dimension lda
 * @details TODO why do we have an r and a c?
 * @param[in] r Description
 * @param[in] c Description
 * @return Description
 */
inline size_t RGF<T>::mf_block_lda(size_t r, size_t c) {
  // return matrix->index_i[c*b_size];
  // two blocks
  if (c < matrix_nt - 1)
    return 2 * matrix_ns + matrix_nd;
  // one block
  if (c < matrix_nt)
    return matrix_ns + matrix_nd;
  // dense block
  else
    return matrix_nd;
}

/************************************************************************************************/

template <class T> inline size_t RGF<T>::mf_dense_block_index(size_t i) {
  if (i < matrix_nt - 1)
    return diag_pos[i * b_size] + 2 * b_size;
  else if (i == matrix_nt - 1)
    return diag_pos[i * b_size] + b_size;
  else
    return diag_pos[i * b_size];
}

/************************************************************************************************/

template <class T> inline size_t RGF<T>::mf_dense_block_offset(size_t i) {
  if (i < matrix_nt - 1)
    return 2 * b_size;
  else if (i == matrix_nt - 1)
    return b_size;
  else
    return 0;
}

/************************************************************************************************/

template <class T> inline size_t RGF<T>::mf_dense_block_lda(size_t i) {
  return mf_block_lda(i, i);
  ;
}

/************************************************************************************************/

template <class T> double RGF<T>::FirstStageFactor() {
  int info;
  /**< TODO IB: index of block*/
  size_t IB;
  /**< NR: Number of the rows for the first block = ns*/
  /**< NM: Number of the rows for the seond block = ns*/
  size_t NR, NM;
  T ONE = f_one();

  cudaStream_t init_stream1, init_stream2, d2h_stream;
  cudaEvent_t h2d_event, computation_event;

  cudaStreamCreate(&init_stream1);
  cudaStreamCreate(&init_stream2);
  cudaStreamCreate(&d2h_stream);
  cudaEventCreate(&h2d_event);
  cudaEventCreate(&computation_event);
  // count floating point operations
  double flops = 0;

  // tpotrf_dev('L', NR, M1_dev ,NR ,&info)
  // rj dpotrf MF[0,0]
  // rj dtrsm RLTN MF[0,0] MF[1,0]

  IB = 0;
  NR = Bmax[0] - Bmin[0];

  // Allocate inital block and if exist fixed effect blocks
  init_supernode(blockR_dev, IB);
  if (matrix_nd > 0) {
    init_supernode(blockDense_dev, matrix_nt);
  }

  // For the first Block /////////////////////////////////////////////////////
  /**< CHOLESKY FACTORISATION A = L*L**T of the first NR=ns*ns diagonal Block
  NAMING: type_positive_definite_triangular_factorization = Cholesky
  factorization*/
  tpotrf_dev('L', NR, blockR_dev, mf_block_lda(IB, IB), &info);
  // printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR, mf_block_index(0, 0),
  // mf_block_lda(0, 0));
  // mf_block_lda(0, 0));

  // TODO N^3/3 FLOPS instead of 0.5
  flops += 1.0 / 3.0 * NR * NR * NR;

  // TODO: shouldn't that be if matrix_nd>1 ?
  // rj dtrsm RLTN MF[IB,IB] MF[IB+1,IB]
  if (matrix_nt > 1) {
    // UPDATE COLUMNS ACCORDING TO CHOLESKY FACTORISATION
    /**< triangular solve matrix: \p A*\p X = \p alpha \p B
     * X, B mxn
     * side: of A L or R
     * uplo: is A upper or lower U or L
     * tran: is A a transponsed or not N or T
     * diag: is A a unit triangular U or N
     * m: number of rows of B
     * n: number columns of B
     * alpha
     * A lda x k
     * k is m if side=L and n if side is R
     * lda: firs dimension of A
     * b: on exit is overwritten by the solution matrix X
     * ldb: first dimension of B at least max(1,m)
     */
    // TODO: calculates C_1 understand this in detail
    ttrsm_dev('R', 'L', 'T', 'N', NR + matrix_nd, NR, ONE, blockR_dev,
              mf_block_lda(0, 0), &blockR_dev[NR], mf_block_lda(1, 0),
              magma_queue);
    flops += 2.0 * (NR + matrix_nd) * NR * NR;
  }

  copy_supernode_to_host(blockR_dev, IB);

  // For the other Blocks /////////////////////////////////////////////////////
  // rj IB = 1; IB < NBlock; IB++
  for (IB = 1; IB < matrix_nt - 1; IB++) {
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    // NOTE:
    // here we swap the blocks A_i<-> TODO
    // In the first round blockM_dev=blockR_dev=C_1
    // and blockR_dev=0
    // CPU call
    // TODO maybecreate many with a for loop
    if(IB % 2 == 1){
      swap_pointers(&blockR_dev, &blockR_dev1);
    } else{
      swap_pointers(&blockR_dev, &blockR_dev2);
    }
    swap_pointers(&blockR_dev, &blockM_dev);

    init_supernode(blockR_dev1, IB, init_stream1);
    init_supernode(blockR_dev2, IB+1, init_stream2);

    /** Note: init_supernode does the following
    // KERNEL
    // fill_dev(M_dev, T(0.0), supernode_size);
    // => Calls a kernel<girdDim, blockdim> to allocate memory for M_dev on dev
    and
    // fills it with zero
    //
    // MEMCPY
    // Copy tiplet pointers to device pointers ia_dev, ja_dev, a_dev
    // memcpy_to_device(&matrix_ia[supernode_fc], ia_dev, (cols + 1) *
    sizeof(size_t));
    // memcpy_to_device(&matrix_ja[supernode_offset], ja_dev, supernode_nnz *
    sizeof(size_t));
    // memcpy_to_device(&matrix_a[supernode_offset], a_dev, supernode_nnz *
    sizeof(T));
    // => they use cudaMemcpyAsync with NULL stream
    //
    // KERNEL
    // set M_dev on dev from ia_dev, ja_dev, a_dev
    // init_supernode_dev(M_dev, ia_dev, ja_dev, a_dev, supernode,
    //                    supernode_nnz, supernode_offset, matrix_ns, matrix_nt,
    matrix_nd);
    // => Calls a kernel<gridDim, blockdim> to allocate memory for M_dev and
    sets the matrix from
    // the the triplets ia_dev, ja_dev, M_dev
    // afterwards calls __syncthreads();
    */
    // UPDATE NEXT DIAGONAL BLOCK
    // rj dgemm NT M[IB,IB-1] M[IB,IB-1]
    // todo rj: 3-last parameter ZERO in PARDISO
    // C_i C_i^T
    /**<
     * S_i=(A_i-C_i S_i B_i)=(A_i-C_i blockM_dev)
     * TODO: why do we have &blockM_dev[NM] and also
     * blockM_dev[NM]*blockM_dev[NM]^T?
     */
    // enqueue event in stream
    if(IB % 2 == 1){
      cudaEventRecord(h2d_event, init_stream1);
    } else{
      cudaEventRecord(h2d_event, init_stream2);
    }
    cudaStreamWaitEvent(NULL, h2d_event, 0);
    /**<
     * C=beta C+alpha*AB
     * ...,alpha,A,...,B,...,beta,C
     */
    tgemm_dev('N', 'T', NR, NR, NM, -ONE, &blockM_dev[NM],
              mf_block_lda(IB, IB - 1), &blockM_dev[NM],
              mf_block_lda(IB, IB - 1), ONE, blockR_dev, mf_block_lda(IB, IB),
              magma_queue);
    flops += 2.0 * NR * NR * NR;

    if (matrix_nd > 0) {
      // Update dense rows of super node IB
      tgemm_dev('N', 'T', matrix_nd, NR, NM, -ONE,
                &blockM_dev[mf_dense_block_offset(IB - 1)],
                mf_dense_block_lda(IB - 1), &blockM_dev[NM],
                mf_block_lda(IB, IB - 1), ONE,
                &blockR_dev[mf_dense_block_offset(IB)], mf_dense_block_lda(IB),
                magma_queue);
      // NM = NR
      flops += 2.0 * matrix_nd * NR * NM;

      // update last diagonal block
      tgemm_dev('N', 'T', matrix_nd, matrix_nd, NM, -ONE,
                &blockM_dev[mf_dense_block_offset(IB - 1)],
                mf_dense_block_lda(IB - 1),
                &blockM_dev[mf_dense_block_offset(IB - 1)],
                mf_dense_block_lda(IB - 1), ONE, blockDense_dev, matrix_nd,
                magma_queue);
      flops += 2.0 * matrix_nd * matrix_nd * NM;
    }
    // printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c:
    // %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB,
    // IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1),
    // mf_block_index(IB, IB), mf_block_lda(IB, IB)); rj dpotrf MF[IB,IB]
    tpotrf_dev('L', NR, blockR_dev, mf_block_lda(IB, IB), &info);
    flops += 0.5 * NR * NR * NR;
    // printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR,
    // mf_block_index(IB, IB), mf_block_lda(IB, IB)); rj dtrsm RLTN
    // MF[IB,IB] MF[IB+1,IB]
    // TODO: check if this is true: NOTE: only the first block matrix is of full
    // size and the others are triangular!
    ttrsm_dev('R', 'L', 'T', 'N', NR + matrix_nd, NR, ONE, blockR_dev,
              mf_block_lda(IB, IB), &blockR_dev[NR], mf_block_lda(IB + 1, IB),
              magma_queue);
    flops += 2.0 * (NR + matrix_nd) * NR * NR;
    // printf("RJ: ttrsm NR: %d, a: %d, lda: %d, b: %d, ldb: %d\n", NR,
    // mf_block_index(IB-1, IB-1), mf_block_lda(IB-1, IB-1),
    // mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1));
    // TODO: can be done during initalizition of the next supernode
    // NOTE: at the moment we just have one stream and just swap the pointers
    // and then initialize If we use two stream and make sure that the pointers
    // are swaped on the device then we can start the initalization while also
    // copying copy_supernode_to_host(blockR_dev, IB, STREAM); -> uses
    // cudaMemcpyAsync already in background with NULL stream enqueue event in
    // stream
    cudaEventRecord(computation_event, NULL);
    cudaStreamWaitEvent(d2h_stream, computation_event, 0);
    copy_supernode_to_host(blockR_dev, IB, d2h_stream);
    /** Note:
    // Calls:
    // memcpy_to_host(&MF[ind], M_dev=blockR_dev, rows * cols * sizeof(T));
    // which in turn calls
    // cudaMemcpyAsync(host_data=&MF[ind], device_data=blockR_dev,
    size_element,cudaMemcpyDeviceToHost, NULL);
    */
  }

  if (matrix_nt > 1) {
    IB = matrix_nt - 1;
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    swap_pointers(&blockR_dev, &blockM_dev);
    init_supernode(blockR_dev, IB);

    // rj dgemm NT M[IB,IB-1] M[IB,IB-1]
    // todo rj: 3-last parameter ZERO in PARDISO
    tgemm_dev('N', 'T', NR, NR, NM, -ONE, &blockM_dev[NM],
              mf_block_lda(IB, IB - 1), &blockM_dev[NM],
              mf_block_lda(IB, IB - 1), ONE, blockR_dev, mf_block_lda(IB, IB),
              magma_queue);
    flops += 2.0 * NR * NR * NR;

    if (matrix_nd > 0) {
      tgemm_dev('N', 'T', matrix_nd, NR, NM, -ONE,
                &blockM_dev[mf_dense_block_offset(IB - 1)],
                mf_dense_block_lda(IB - 1), &blockM_dev[NM],
                mf_block_lda(IB, IB - 1), ONE,
                &blockR_dev[mf_dense_block_offset(IB)], mf_dense_block_lda(IB),
                magma_queue);
      flops += 2.0 * matrix_nd * NR * NM;

      tgemm_dev('N', 'T', matrix_nd, matrix_nd, NM, -ONE,
                &blockM_dev[mf_dense_block_offset(IB - 1)],
                mf_dense_block_lda(IB - 1),
                &blockM_dev[mf_dense_block_offset(IB - 1)],
                mf_dense_block_lda(IB - 1), ONE, blockDense_dev, matrix_nd,
                magma_queue);
      flops += 2.0 * matrix_nd * matrix_nd * NM;
    }
    // printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c:
    // %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB,
    // IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1),
    // mf_block_index(IB, IB), mf_block_lda(IB, IB)); rj dpotrf MF[IB,IB]
    tpotrf_dev('L', NR, blockR_dev, mf_block_lda(IB, IB), &info);
    flops += 0.5 * NR * NR * NR;

    // printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR,
    // mf_block_index(IB, IB), mf_block_lda(IB, IB)); rj dtrsm RLTN
    // MF[IB,IB] MF[IB+1,IB]
    if (matrix_nd > 0) {
      ttrsm_dev('R', 'L', 'T', 'N', matrix_nd, NR, ONE, blockR_dev,
                mf_block_lda(IB, IB), &blockR_dev[mf_dense_block_offset(IB)],
                mf_dense_block_lda(IB), magma_queue);
      flops += 2.0 * matrix_nd * NR * NR;
    }
    // printf("RJ: ttrsm NR: %d, a: %d, lda: %d, b: %d, ldb: %d\n", NR,
    // mf_block_index(IB-1, IB-1), mf_block_lda(IB-1, IB-1),
    // mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1));

    copy_supernode_to_host(blockR_dev, IB);
    // TODO for copying non-zero blocks only
    // CudaMemcpyAsync blockR_dev->tmp_buffer
    // returns
    // blockRdev
    // CudaMemcpyAsync is copying tmp_buffer
  }

  if (matrix_nd > 0) {
    IB = NBlock - 1;
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    // rj dgemm NT M[IB,IB-1] M[IB,IB-1]
    // rj NM is the same for all blocks 0..nt-1
    tgemm_dev(
        'N', 'T', matrix_nd, matrix_nd, NM, -ONE,
        &blockR_dev[mf_dense_block_offset(IB - 1)], mf_dense_block_lda(IB - 1),
        &blockR_dev[mf_dense_block_offset(IB - 1)], mf_dense_block_lda(IB - 1),
        ONE, blockDense_dev, matrix_nd, magma_queue);
    flops += 2.0 * matrix_nd * matrix_nd * NM;

    // printf("RJ: tgemm NR: %d, NM: %d, a: %d, lda: %d, b: %d, ldb: %d, c:
    // %d, ldc: %d\n", NR, NM, mf_block_index(IB, IB-1), mf_block_lda(IB,
    // IB-1), mf_block_index(IB, IB-1), mf_block_lda(IB, IB-1),
    // mf_block_index(IB, IB), mf_block_lda(IB, IB)); rj dpotrf MF[IB,IB]
    tpotrf_dev('L', NR, blockDense_dev, mf_block_lda(IB, IB), &info);
    flops += 0.5 * NR * NR * NR;

    // printf("RJ: tpotrf NR: %d, a: %d, lda: %d\n\n", NR,
    // mf_block_index(IB, IB), mf_block_lda(IB, IB));

    copy_supernode_to_host(blockDense_dev, IB);
  }

  cudaStreamDestroy(init_stream1);
  cudaStreamDestroy(init_stream2);
  cudaStreamDestroy(d2h_stream);
  cudaEventDestroy(h2d_event);
  cudaEventDestroy(computation_event);

  return flops;
}

/************************************************************************************************/
/**
 * @brief solving of the linear system Qx=b
 * @details Description
 * \todo TODO: slow due to many zeros
 * Question: can we learn where we have the zeros?
 * @param[in] nrhs Description
 */
template <class T> double RGF<T>::SecondStageSolve(size_t nrhs) {
  size_t IB;
  size_t NR, NM; // NP;
  T ONE = f_one();

  // Forward pass
  NR = Bmax[0] - Bmin[0];

  // counting FLOPS
  double flops = 0;

  // rj dtrsm LLNN MF[0,0] rhs[0]
  // solve for first diagonal block
  c_ttrsm('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF[mf_block_index(0, 0)],
          mf_block_lda(0, 0), &rhs[Bmin[0]], matrix_size);
  flops += 2.0 * NR * NR * nrhs;

  // rj IB = 1; IB < NBlock; IB++
  for (IB = 1; IB < matrix_nt; IB++) {
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    // rj dgemm NN M[IB,IB-1] rhs[IB-1] rhs[IB]
    c_tgemm('N', 'N', NR, nrhs, NM, -ONE, &MF[mf_block_index(IB, IB - 1)],
            mf_block_lda(IB, IB - 1), &rhs[Bmin[IB - 1]], matrix_size, ONE,
            &rhs[Bmin[IB]], matrix_size);
    flops += 2.0 * NR * nrhs * NM;
    // rj dtrsm LLNN MF[IB,IB] rhs[IB]
    c_ttrsm('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF[mf_block_index(IB, IB)],
            mf_block_lda(IB, IB), &rhs[Bmin[IB]], matrix_size);
    flops += 2.0 * NR * NR * nrhs;
  }

  if (matrix_nd > 0) {
    IB = NBlock - 1;
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    // rj dgemm NN M[IB,IB-1] rhs[IB-1] rhs[IB]
    for (size_t i = 0; i < NBlock - 1; i++) {
      c_tgemm('N', 'N', NR, nrhs, NM, -ONE, &MF[mf_dense_block_index(i)],
              mf_dense_block_lda(i), &rhs[Bmin[i]], matrix_size, ONE,
              &rhs[Bmin[IB]], matrix_size);
      flops += 2.0 * NR * nrhs * NM;
    }
    // rj dtrsm LLNN MF[IB,IB] rhs[IB]
    c_ttrsm('L', 'L', 'N', 'N', NR, nrhs, ONE, &MF[mf_block_index(IB, IB)],
            mf_block_lda(IB, IB), &rhs[Bmin[IB]], matrix_size);
    flops += 2.0 * NR * NR * nrhs;
  }

  // Backward pass
  if (matrix_nd > 0) {
    IB = NBlock - 1;
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    // rj dtrsm LLTN MF[NBlock-1,NBlock-1] rhs[NBlock-1]
    c_ttrsm('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF[mf_block_index(IB, IB)],
            mf_block_lda(IB, IB), &rhs[Bmin[IB]], matrix_size);
    flops += 2.0 * NR * NR * nrhs;
    // rj dgemm TN M[IB,IB-1] rhs[IB] rhs[IB-1]
    for (size_t i = 0; i < NBlock - 1; i++) {
      c_tgemm('T', 'N', NM, nrhs, NR, -ONE, &MF[mf_dense_block_index(i)],
              mf_dense_block_lda(i), &rhs[Bmin[IB]], matrix_size, ONE,
              &rhs[Bmin[i]], matrix_size);
      flops += 2.0 * NM * nrhs * NR;
    }
  }

  for (IB = matrix_nt - 1; IB > 0; IB--) {
    NR = Bmax[IB] - Bmin[IB];
    NM = Bmax[IB - 1] - Bmin[IB - 1];

    // rj dtrsm LLTN MF[IB,IB] rhs[IB]
    c_ttrsm('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF[mf_block_index(IB, IB)],
            mf_block_lda(IB, IB), &rhs[Bmin[IB]], matrix_size);
    flops += 2.0 * NR * NR * nrhs;
    // rj dgemm TN M[IB,IB-1] rhs[IB] rhs[IB-1]
    c_tgemm('T', 'N', NM, nrhs, NR, -ONE, &MF[mf_block_index(IB, IB - 1)],
            mf_block_lda(IB, IB - 1), &rhs[Bmin[IB]], matrix_size, ONE,
            &rhs[Bmin[IB - 1]], matrix_size);
    flops += 2.0 * NM * nrhs * NR;
  }

  IB = 0;
  NR = Bmax[IB] - Bmin[IB];

  // rj dtrsm LLTN MF[NBlock-1,NBlock-1] rhs[NBlock-1]
  c_ttrsm('L', 'L', 'T', 'N', NR, nrhs, ONE, &MF[mf_block_index(IB, IB)],
          mf_block_lda(IB, IB), &rhs[Bmin[IB]], matrix_size);
  flops += 2.0 * NR * NR * nrhs;

  return flops;
}

/************************************************************************************************/
/**
 * @brief Calculates G_i i=n+1,...,1
 * @details Description
 * @param[inout] tmp1_dev Description
 * @param[out] tmp2_dev Description
 */
template <class T> double RGF<T>::ThirdStageRGF(T *tmp1_dev, T *tmp2_dev) {
  int info;
  size_t IB;
  size_t NR, NM, NP, ND;
  cudaStream_t init_stream1, init_stream2;
  cudaStream_t h2d_stream;
  cudaEvent_t h2d_event, computation_event;
  T ONE = f_one();
  T ZERO = f_zero();

  ND = matrix_nd;

  double flops = 0;
  cudaStreamCreate(&init_stream1);
  cudaStreamCreate(&init_stream2);
  cudaStreamCreate(&h2d_stream);
  cudaEventCreate(&h2d_event);
  cudaEventCreate(&computation_event);

  if (matrix_nd > 0) {
    // dense block
    IB = NBlock - 1;
    NR = Bmax[IB] - Bmin[IB];

    tlacpy_dev('L', NR, NR, blockDense_dev, mf_block_lda(IB, IB), tmp1_dev, NR,
               magma_queue);
    tril_dev(tmp1_dev, NR, NR);
    ttrtri_dev('L', 'N', NR, tmp1_dev, NR, &info);
    flops += 0.5 * NR * NR * NR;
    tgemm_dev('T', 'N', NR, NR, NR, ONE, tmp1_dev, NR, tmp1_dev, NR, ZERO,
              blockDense_dev, mf_block_lda(IB, IB), magma_queue);
    flops += 2.0 * NR * NR * NR;

    // last non-dense block
    IB = NBlock - 2;
    NR = Bmax[IB] - Bmin[IB];
    NP = Bmax[IB + 1] - Bmin[IB + 1];

    tril_dev(blockR_dev, mf_block_lda(IB, IB), NR);
    ttrtri_dev('L', 'N', NR, blockR_dev, mf_block_lda(IB, IB), &info);
    flops += 0.5 * NR * NR * NR;
    init_eye_on_dev(tmp1_dev, NR, 0);

    tgemm_dev('T', 'N', NR, NP, NP, ONE, &blockR_dev[mf_dense_block_offset(IB)],
              mf_dense_block_lda(IB), blockDense_dev,
              mf_block_lda(IB + 1, IB + 1), ZERO, tmp2_dev, NR, magma_queue);
    flops += 2.0 * NR * NP * NP;
    tgemm_dev('N', 'N', NR, NR, NP, ONE, tmp2_dev, NR,
              &blockR_dev[mf_dense_block_offset(IB)], mf_dense_block_lda(IB),
              ONE, tmp1_dev, NR, magma_queue);
    flops += 2.0 * NR * NR * NP;
    tgemm_dev('T', 'N', NR, NR, NR, ONE, blockR_dev, mf_block_lda(IB, IB),
              tmp1_dev, NR, ZERO, tmp2_dev, NR, magma_queue);
    flops += 2.0 * NR * NR * NR;
    tgemm_dev('N', 'N', NR, NR, NR, ONE, tmp2_dev, NR, blockR_dev,
              mf_block_lda(IB, IB), ZERO, tmp1_dev, NR, magma_queue);
    flops += 2.0 * NR * NR * NR;
    tlacpy_dev('F', NR, NR, tmp1_dev, NR, blockR_dev, mf_block_lda(IB, IB),
               magma_queue);

    copy_supernode_diag(blockDense_dev, NBlock - 1);
  } else {
    // last non-dense block
    IB = NBlock - 1;
    NR = Bmax[IB] - Bmin[IB];

    tlacpy_dev('L', NR, NR, blockR_dev, mf_block_lda(IB, IB), tmp1_dev, NR,
               magma_queue);
    tril_dev(tmp1_dev, NR, NR);
    ttrtri_dev('L', 'N', NR, tmp1_dev, NR, &info);
    flops += 0.5 * NR * NR * NR;
    tgemm_dev('T', 'N', NR, NR, NR, ONE, tmp1_dev, NR, tmp1_dev, NR, ZERO,
              blockR_dev, mf_block_lda(IB, IB), magma_queue);
  }

  // second-last non-dense block .. 0-block
  for (int IBi = matrix_nt - 2; IBi > -1; IBi--) {
    IB = IBi;
    NR = Bmax[IB] - Bmin[IB];
    NP = Bmax[IB + 1] - Bmin[IB + 1];

    // TODO
    swap_pointers(&blockR_dev, &blockP_dev);

    // enqueue event in stream
    cudaEventRecord(computation_event, NULL);
    cudaStreamWaitEvent(h2d_stream, computation_event, 0);

    copy_supernode_diag(blockP_dev, IB + 1);
    copy_supernode_to_device(blockR_dev, IB, h2d_stream);

    // enqueue event in stream
    cudaEventRecord(h2d_event, h2d_stream);
    cudaStreamWaitEvent(NULL, h2d_event, 0);

    tril_dev(blockR_dev, mf_block_lda(IB, IB), NR);
    ttrtri_dev('L', 'N', NR, blockR_dev, mf_block_lda(IB, IB), &info);
    flops += 0.5 * NR * NR * NR;
    init_eye_on_dev(tmp1_dev, NR, 0);

    tgemm_dev('T', 'N', NR, NP, NP, ONE, &blockR_dev[NR],
              mf_block_lda(IB + 1, IB), blockP_dev,
              mf_block_lda(IB + 1, IB + 1), ZERO, tmp2_dev, NR, magma_queue);
    flops += 2.0 * NR * NP * NP;
    tgemm_dev('N', 'N', NR, NR, NP, ONE, tmp2_dev, NR, &blockR_dev[NR],
              mf_block_lda(IB + 1, IB), ONE, tmp1_dev, NR, magma_queue);
    flops += 2.0 * NR * NR * NP;

    if (matrix_nd > 0) {
      tgemm_dev('T', 'N', NR, ND, ND, ONE,
                &blockR_dev[mf_dense_block_offset(IB)], mf_dense_block_lda(IB),
                blockDense_dev, mf_block_lda(NBlock - 1, NBlock - 1), ZERO,
                tmp2_dev, NR, magma_queue);
      flops += 2.0 * NR * ND * ND;
      tgemm_dev('N', 'N', NR, NR, ND, ONE, tmp2_dev, NR,
                &blockR_dev[mf_dense_block_offset(IB)], mf_dense_block_lda(IB),
                ONE, tmp1_dev, NR, magma_queue);
      flops += 2.0 * NR * NR * ND;
    }

    tgemm_dev('T', 'N', NR, NR, NR, ONE, blockR_dev, mf_block_lda(IB, IB),
              tmp1_dev, NR, ZERO, tmp2_dev, NR, magma_queue);
    flops += 2.0 * NR * NR * NR;
    tgemm_dev('N', 'N', NR, NR, NR, ONE, tmp2_dev, NR, blockR_dev,
              mf_block_lda(IB, IB), ZERO, tmp1_dev, NR, magma_queue);
    flops += 2.0 * NR * NR * NR;
    tlacpy_dev('F', NR, NR, tmp1_dev, NR, blockR_dev, mf_block_lda(IB, IB),
               magma_queue);
  }

  copy_supernode_diag(blockR_dev, 0);

  return flops;
}

/************************************************************************************************/
/**
 * @brief Blocked Cholesky Factorization
 * @details Compute the cholesky factor for all the blocks uning block Cholesky
 * factorization is the second step after the constuctor.
 */
template <class T> double RGF<T>::factorize() {
  ///
  /// finds maximum over all blocks and set it to b_size
  /// for us it should be either ns*ns or nd*ns if nd > ns
  create_blocks();

  // Data allocation
  if (!MF_dev_allocated) {
    MF = new T[matrix_nnz];
    /**< max_supernode_nnz: size of supernode blocks
     * for nt=0 we have just one single block
     * for nt>0 we have a big supernode block
     * NOTE: a supernode block [A B fixed]^T
     */
    size_t max_supernode_nnz = matrix_nt > 1
                                   ? matrix_ns * (2 * matrix_ns + matrix_nd)
                                   : matrix_ns * (matrix_ns + matrix_nd);
    /**< dense_supernode_nnz: size of dense random effects blocks if we have
     * fixed effects*/
    size_t dense_supernode_nnz = matrix_nd > 0 ? matrix_nd * matrix_nd : 0;
    allocate_data_on_device((void **)&blockR_dev,
                            max_supernode_nnz * sizeof(T));
    allocate_data_on_device((void **)&blockR_dev1,
                            max_supernode_nnz * sizeof(T));
    allocate_data_on_device((void **)&blockR_dev2,
                            max_supernode_nnz * sizeof(T));
    allocate_data_on_device((void **)&blockM_dev,
                            max_supernode_nnz * sizeof(T));
    allocate_data_on_device((void **)&blockP_dev,
                            max_supernode_nnz * sizeof(T));
    allocate_data_on_device((void **)&blockDense_dev,
                            dense_supernode_nnz * sizeof(T));
    MF_dev_allocated = true;
  }

  // size_t nnz = matrix_ia[matrix_size];
  /**< max_rows: number of supernode blocks
   * for nt=0 we have just one single block
   * for nt>0 we have a big supernode block
   */
  size_t max_rows =
      matrix_nt > 1 ? 2 * matrix_ns + matrix_nd : matrix_ns + matrix_nd;
  size_t max_cols = fmax(matrix_ns, matrix_nd);

  // Temp data allocation for triplet format
  allocate_data_on_device((void **)&ia_dev, (max_cols + 1) * sizeof(size_t));
  allocate_data_on_device((void **)&ja_dev,
                          max_rows * max_cols * sizeof(size_t));
  allocate_data_on_device((void **)&a_dev, max_rows * max_cols * sizeof(T));

  // Computation
  double flops = FirstStageFactor();

  // Temp data deallocation
  deallocate_data_on_dev(ia_dev, (max_cols + 1) * sizeof(size_t));
  deallocate_data_on_dev(ja_dev, max_rows * max_cols * sizeof(size_t));
  deallocate_data_on_dev(a_dev, max_rows * max_cols * sizeof(T));

  factorization_completed = true;

  return flops;
}

/************************************************************************************************/

template <class T> double RGF<T>::solve(T *b, size_t nrhs) {
  double flops = solve(b, b, nrhs);

  return flops;
}

/************************************************************************************************/
/**
 * @brief Summary
 * @details Description
 * @param[inout] x Description
 * @param[out] b Description
 * @param[in] nrhs Description
 */
template <class T> double RGF<T>::solve(T *x, T *b, size_t nrhs) {
  if (!factorization_completed)
    factorize();

  // Data allocation
  rhs = new T[matrix_size * nrhs];

  // Copy data
  memcpy(rhs, b, matrix_size * nrhs * sizeof(T));

  // Computation
  double flops = SecondStageSolve(nrhs);

  // Copy data
  memcpy(x, rhs, matrix_size * nrhs * sizeof(T));

  // Data deallocation
  delete[] rhs;

  return flops;
}

/************************************************************************************************/

template <class T> double RGF<T>::RGFdiag(T *diag) {
  if (!factorization_completed)
    factorize();

  factorization_completed = false;

  // Data allocation
  T *tmp1_dev;
  T *tmp2_dev;
  allocate_data_on_device((void **)&diag_dev, matrix_size * sizeof(T));
  allocate_data_on_device((void **)&diag_pos_dev, matrix_size * sizeof(size_t));
  allocate_data_on_device((void **)&tmp1_dev, b_size * b_size * sizeof(T));
  allocate_data_on_device((void **)&tmp2_dev, b_size * b_size * sizeof(T));

  // Copy diag_pos to device
  memcpy_to_device(diag_pos, diag_pos_dev, matrix_size * sizeof(size_t));

  // Computation
  double flops = ThirdStageRGF(tmp1_dev, tmp2_dev);

  // Copy data to host
  memcpy_to_host(diag, diag_dev, matrix_size * sizeof(T));

  // Data deallocation
  deallocate_data_on_dev(diag_dev, matrix_size * sizeof(T));
  deallocate_data_on_dev(diag_pos_dev, matrix_size * sizeof(size_t));
  deallocate_data_on_dev(tmp1_dev, b_size * b_size * sizeof(T));
  deallocate_data_on_dev(tmp2_dev, b_size * b_size * sizeof(T));

  return flops;
}

/************************************************************************************************/

template <class T> T RGF<T>::logDet() {
  if (!factorization_completed)
    factorize();

  // Computation
  T det = T(0.0);
  indexed_log_sum(MF, diag_pos, matrix_size, &det);

  return 2.0 * det;
}

/************************************************************************************************/

template <class T> double RGF<T>::residualNorm(T *x, T *b) {
  T *r = new T[matrix_size];

  memcpy(r, b, matrix_size * sizeof(T));

  for (size_t ic = 0; ic < matrix_size; ic++) {
    for (size_t i = matrix_ia[ic]; i < matrix_ia[ic + 1]; i++) {
      size_t ir = matrix_ja[i];

      r[ir] -= matrix_a[i] * x[ic];
      if (ir != ic)
        r[ic] -= matrix_a[i] * x[ir];
    }
  }

  double res = c_dtnrm2(matrix_size, r, 1);

  delete[] r;

  return res;
}

/************************************************************************************************/

template <class T> double RGF<T>::residualNormNormalized(T *x, T *b) {
  return residualNorm(x, b) / c_dtnrm2(matrix_size, b, 1);
}

/************************************************************************************************/

/**
 * @brief Calculates the largest block size
 * @details finds maximum over all blocks and sets \p b_size to it.
 * Note: for our application the b_size is const apart from the fixed effect
 * blocks.
 */
template <class T> void RGF<T>::create_blocks() {
  size_t IB;

  b_size = 0;

  for (IB = 0; IB < NBlock; IB++) {

    if (Bmax[IB] - Bmin[IB] > b_size) {
      b_size = Bmax[IB] - Bmin[IB];
    }
  }
}

/************************************************************************************************/

template <class T>

/**
 * @brief Select a supernod that we want to initialize
 * Supernode = A block together with random effects
 * @details initializes M_dev supernode from a
 * @param[inout] M_dev supernode on device to be initalized
 * @param[in] supernode number of supernode to be initalize
 */
inline void RGF<T>::init_supernode(T *M_dev, size_t supernode) {
  size_t supernode_fc = supernode * matrix_ns;
  size_t supernode_lc = supernode < matrix_nt
                            ? (supernode + 1) * matrix_ns
                            : matrix_ns * matrix_nt + matrix_nd;
  size_t supernode_nnz = matrix_ia[supernode_lc] - matrix_ia[supernode_fc];
  size_t supernode_offset =
      matrix_ia[supernode_fc]; /**< start of the first supernode*/
  size_t rows = mf_block_lda(supernode, supernode); /**< rows of the supernode*/
  size_t cols = supernode_lc - supernode_fc;        /**< cols of the supernode*/
  size_t supernode_size = rows * cols;

  fill_dev(M_dev, T(0.0), supernode_size);
  // Copy tiplet pointers to device
  memcpy_to_device(&matrix_ia[supernode_fc], ia_dev,
                   (cols + 1) * sizeof(size_t));
  memcpy_to_device(&matrix_ja[supernode_offset], ja_dev,
                   supernode_nnz * sizeof(size_t));
  memcpy_to_device(&matrix_a[supernode_offset], a_dev,
                   supernode_nnz * sizeof(T));
  // use triplet pointers on device to set M_dev
  init_supernode_dev(M_dev, ia_dev, ja_dev, a_dev, supernode, supernode_nnz,
                     supernode_offset, matrix_ns, matrix_nt, matrix_nd);
}

/************************************************************************************************/

template <class T>

/**
 * @brief Select a supernod that we want to initialize
 * @details initializes M_dev supernode from a
 * @param[inout] M_dev supernode on device to be initalized
 * @param[in] supernode number of supernode to be initalize
 */
inline void RGF<T>::init_supernode(T *M_dev, size_t supernode,
                                   cudaStream_t stream) {
  size_t supernode_fc = supernode * matrix_ns;
  size_t supernode_lc = supernode < matrix_nt
                            ? (supernode + 1) * matrix_ns
                            : matrix_ns * matrix_nt + matrix_nd;
  size_t supernode_nnz = matrix_ia[supernode_lc] - matrix_ia[supernode_fc];
  size_t supernode_offset = matrix_ia[supernode_fc];
  size_t rows = mf_block_lda(supernode, supernode);
  size_t cols = supernode_lc - supernode_fc;
  size_t supernode_size = rows * cols;

  fill_dev(M_dev, T(0.0), supernode_size, stream);
  // Copy tiplet pointers to device
  memcpy_to_device(&matrix_ia[supernode_fc], ia_dev,
                   (cols + 1) * sizeof(size_t), stream);
  memcpy_to_device(&matrix_ja[supernode_offset], ja_dev,
                   supernode_nnz * sizeof(size_t), stream);
  memcpy_to_device(&matrix_a[supernode_offset], a_dev,
                   supernode_nnz * sizeof(T), stream);
  // use triplet pointers on device to set M_dev
  init_supernode_dev(M_dev, ia_dev, ja_dev, a_dev, supernode, supernode_nnz,
                     supernode_offset, matrix_ns, matrix_nt, matrix_nd, stream);
}

/************************************************************************************************/

template <class T>
inline void RGF<T>::copy_supernode_to_host(T *M_dev, size_t supernode,
                                           cudaStream_t stream) {
  size_t supernode_fc = supernode * matrix_ns;
  size_t supernode_lc = supernode < matrix_nt
                            ? (supernode + 1) * matrix_ns
                            : matrix_ns * matrix_nt + matrix_nd;
  size_t ind = mf_block_index(supernode, supernode);
  size_t rows = mf_block_lda(supernode, supernode);
  size_t cols = supernode_lc - supernode_fc;

  memcpy_to_host(&MF[ind], M_dev, rows * cols * sizeof(T), stream);
}

/************************************************************************************************/

template <class T>
inline void RGF<T>::copy_supernode_to_device(T *M_dev, size_t supernode,
                                             cudaStream_t stream) {
  size_t supernode_fc = supernode * matrix_ns;
  size_t supernode_lc = supernode < matrix_nt
                            ? (supernode + 1) * matrix_ns
                            : matrix_ns * matrix_nt + matrix_nd;
  size_t ind = mf_block_index(supernode, supernode);
  size_t rows = mf_block_lda(supernode, supernode);
  size_t cols = supernode_lc - supernode_fc;

  memcpy_to_device(&MF[ind], M_dev, rows * cols * sizeof(T), stream);
}

/************************************************************************************************/

template <class T>
inline void RGF<T>::copy_supernode_diag(T *src, size_t supernode,
                                        cudaStream_t stream) {
  size_t supernode_fc = supernode * matrix_ns;
  size_t supernode_lc = supernode < matrix_nt
                            ? (supernode + 1) * matrix_ns
                            : matrix_ns * matrix_nt + matrix_nd;
  size_t offset = mf_block_index(supernode, supernode);
  // size_t rows = mf_block_lda(supernode, supernode);
  size_t cols = supernode_lc - supernode_fc;

  indexed_copy_offset_dev(src, &diag_dev[supernode_fc],
                          &diag_pos_dev[supernode_fc], cols, offset, stream);
}

/************************************************************************************************/

template <class T> inline void RGF<T>::swap_pointers(T **ptr1, T **ptr2) {
  T *tmp = *ptr1;
  *ptr1 = *ptr2;
  *ptr2 = tmp;
}

/************************************************************************************************/

template <> CPX RGF<CPX>::f_one() { return CPX(1.0, 0.0); }

template <> double RGF<double>::f_one() { return 1.0; }

/************************************************************************************************/

template <> CPX RGF<CPX>::f_zero() { return CPX(0.0, 0.0); }

template <> double RGF<double>::f_zero() { return 0.0; }

/************************************************************************************************/

#endif
