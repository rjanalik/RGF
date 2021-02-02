#ifndef __RGF
#define __RGF

#include <string.h>

#include "cublas_v2.h"
#include "CSR.H"
#include "Types.H"
#include "Utilities.H"
#include "CWC_utility.H"

class RGF{
	
public:

    RGF(TCSR<CPX>*);

    ~RGF();
    void solve_equation(CPX*,int*,int*,int);
				
private:

    TCSR<CPX>* matrix;
    int findx;
    int b_size;

    cudaStream_t stream_c,stream_m;
    void *cublas_handle;
    CPX *M_dev;
    CPX *M1_dev;
    CPX *M2_dev;
    CPX *M3_dev;
    CPX *nnz_dev;
    int *edge_i_dev;
    int *index_j_dev;

    void extract_diag(int*,int*,CPX*,CPX*,int*,int*,int,cudaStream_t);
    void extract_not_diag(int*,int*,CPX*,CPX*,int*,int*,int,int,cudaStream_t);
    void FirstStage(int*,int*,int);
    void SecondStage(CPX*,int*,int*,int);
    void create_blocks(int*,int*,int);
    void write_matrix(const char*,CPX*,int,int);
};

#endif
