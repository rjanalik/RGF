#include "RGF.H"

/************************************************************************************************/

RGF::RGF(TCSR<CPX>* mat)
{
    matrix = mat;
    findx  = 0;
}

/************************************************************************************************/

RGF::~RGF()
{
}

/************************************************************************************************/

void RGF::extract_diag(int *edge_i,int *index_j,CPX *nnz,CPX *D,int *Bmin,int *Bmax,\
		       int index,cudaStream_t stream)
{
    int NR;
    int imin,imax;

    NR   = Bmax[index]-Bmin[index];
    imin = Bmin[index];
    imax = Bmax[index];

    z_init_var_on_dev(D,NR*NR,stream);

    z_extract_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,matrix->first_row,matrix->findx,\
			  stream);
}

/************************************************************************************************/

void RGF::extract_not_diag(int *edge_i,int *index_j,CPX *nnz,CPX *D,int *Bmin,int *Bmax,\
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

    z_init_var_on_dev(D,NR*NC,stream);

    z_extract_not_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,jmin,side,matrix->first_row,\
			      matrix->findx,stream);
}

/************************************************************************************************/

void RGF::FirstStage(int *b_min,int *b_max,int NBlock)
{
    int info;
    int IB;
    int NR,NM,NP;
    int NS1,NS2;
    CPX ONE   = CPX(1.0,0.0);
    CPX ZERO  = CPX(0.0,0.0);
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
	zaxpy_on_dev(cublas_handle,NR*NR,-ONE,&M_dev[NS1*b_size],1,M1_dev,1);
       
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

void RGF::SecondStage(CPX *GR,int *b_min,int *b_max,int NBlock)
{

    MPI_Status status;
    int info;
    int IB;
    int NR,NM,NP;
    int NS1,NS2;
    int N1,NN;
    CPX ONE      = CPX(1.0,0.0);
    CPX ZERO     = CPX(0.0,0.0);
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
    cudaMemcpy(&GR[b_min[IB]*b_size],&M_dev[NS1*b_size],NR*NR*sizeof(CPX),\
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
	zaxpy_on_dev(cublas_handle,NR*NR,ONE,M2_dev,1,&M_dev[NS1*b_size],1);
	    
	//copy data to host
	cudaMemcpy(&GR[b_min[IB]*b_size],&M_dev[NS1*b_size],NR*NR*sizeof(CPX),\
		   cudaMemcpyDeviceToHost);
    }

    delete[] ipiv;
}

/************************************************************************************************/

void RGF::solve_equation(CPX *GR,int *Bmin,int *Bmax,int NBlock)
{
    cudaStreamCreate(&stream_c);
    cudaStreamCreate(&stream_m);

    magma_init();
    cublas_init(&cublas_handle);

    create_blocks(Bmin,Bmax,NBlock);
    
    //Data allocation
    allocate_data_on_device((void**)&M_dev,b_size*matrix->size*sizeof(CPX));
    allocate_data_on_device((void**)&M1_dev,b_size*b_size*sizeof(CPX));
    allocate_data_on_device((void**)&M2_dev,b_size*b_size*sizeof(CPX));
    allocate_data_on_device((void**)&M3_dev,b_size*b_size*sizeof(CPX));

    allocate_data_on_device((void**)&edge_i_dev,(matrix->size+1)*sizeof(int));
    allocate_data_on_device((void**)&index_j_dev,matrix->n_nonzeros*sizeof(int));
    allocate_data_on_device((void**)&nnz_dev,matrix->n_nonzeros*sizeof(CPX));

    memcpy_to_device(matrix->edge_i,edge_i_dev,(matrix->size+1)*sizeof(int));
    memcpy_to_device(matrix->index_j,index_j_dev,matrix->n_nonzeros*sizeof(int));
    memcpy_to_device(matrix->nnz,nnz_dev,matrix->n_nonzeros*sizeof(CPX));

    //Computation
    init_var(GR,Bmax[NBlock-1]*b_size);
    
    FirstStage(Bmin,Bmax,NBlock);
    SecondStage(GR,Bmin,Bmax,NBlock);

    //Data deallocation
    deallocate_data_on_dev(edge_i_dev,(matrix->size+1)*sizeof(int));
    deallocate_data_on_dev(index_j_dev,matrix->n_nonzeros*sizeof(int));
    deallocate_data_on_dev(nnz_dev,matrix->n_nonzeros*sizeof(CPX));

    deallocate_data_on_dev(M_dev,b_size*matrix->size*sizeof(CPX));
    deallocate_data_on_dev(M1_dev,b_size*b_size*sizeof(CPX));
    deallocate_data_on_dev(M2_dev,b_size*b_size*sizeof(CPX));
    deallocate_data_on_dev(M3_dev,b_size*b_size*sizeof(CPX));

    magma_finalize();
    cublas_finalize(cublas_handle);

    cudaStreamDestroy(stream_c);
    cudaStreamDestroy(stream_m);
}

/************************************************************************************************/

void RGF::create_blocks(int *Bmin,int *Bmax,int NBlock)
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

void RGF::write_matrix(const char *filename,CPX *matrix,int NR,int NC)
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
