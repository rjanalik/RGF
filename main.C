#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>

// include EIGEN
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

// include read/write functions
#include <armadillo>
//#include "read_write_functions.cpp"

// generate data
#include "generate_matrices/generate_poisson_mat.cpp"

// RGF solver
#include "RGF.H"

using namespace std;

using Eigen::MatrixXd;
typedef Eigen::VectorXd Vect;
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif

int main(int argc, char *argv[])
{
   if(argc != 1 + 3){
      std::cout << "wrong number of input parameters. " << std::endl;

      std::cerr << "RGF Call : ns nb noRhs" << std::endl;

      std::cerr << "[integer:ns]                number of spatial grid points " << std::endl;
      //std::cerr << "[integer:nt]                number of temporal grid points " << std::endl;
      std::cerr << "[integer:nb]                number of fixed effects" << std::endl;
      std::cerr << "[integer:noRhs]             number of right hand sides" << std::endl;

      exit(1);
   }

   size_t i; // iteration variable
   //std::cout << "reading in example. " << std::endl;

   size_t ns = atoi(argv[1]);
   size_t nt = ns;
   //size_t nt = atoi(argv[2]);
   size_t nb = atoi(argv[2]);
   size_t nrhs = atoi(argv[3]); 
   
   size_t n  = ns*nt + nb;

   std::cout << "n = " << n << ", ns = " << ns << ", nt = " << nt << ", nb = " << nb << ", nrhs = " << nrhs << "\n" << std::endl;

   if(ns < nb){
      std::cerr << "INVALID PARAMETER CONFIGURATION! ns >= nb required!" << std::endl;
      exit(1);
   }

   std::cout << "========= Computing meomory requirements ===========" << std::endl;

   // compute memory requirements to give warning when too big for GPU
   // we have  2*(nt-1)*ns*ns + ns*ns + n*nb non-zeros
   // for rhs we have nrhs*n non-zeros and then if we also store b we have this two times
   // hence : (2*(nt-1)*ns*ns + ns*ns + n*nb + 2*nrhs*n) * 8 / 10^9
   double mem_gb = (2*(nt-1)*ns*ns + ns*ns + n*nb + 2*nrhs*n) * sizeof(T) / pow(10.0,9.0);
   printf("Memory Usage = %f GB\n\n", mem_gb);


   // set nt = 1 if ns > 0 & nt = 0
   if(ns > 0 && nt == 0){
      nt = 1;
   } 

   // also save as string
   std::string ns_s = std::to_string(ns);
   std::string nt_s = std::to_string(nt);
   std::string nb_s = std::to_string(nb);
   //std::string no_s = std::to_string(no); 

   /* ---------------- generate matrix & rhs ---------------- */
   SpMat Q(n,n);
   generate_adapt_poisson_mat(Q, ns, nb);
   //std::cout << "Q:\n" << MatrixXd(Q) << std::endl;

   MatrixXi B(n, nrhs);
   generate_int_rhs(B, n, nrhs);
   //std::cout << "B = " << B.transpose() << std::endl;

   // =========================================================================== //
   //std::cout << "Converting Eigen Matrices to CSR format. " << std::endl;

   // only take lower triangular part of A
   SpMat Q_lower = Q.triangularView<Lower>(); 
   size_t nnz = Q_lower.nonZeros();

   double t0, t1, t2, t3; // to store times
   size_t* ia; 
   size_t* ja;
   T* a; 
   T *b;
   T *x;
   T *invDiag;

   b  = new T[n*nrhs];
   x  = new T[n*nrhs];
   invDiag = new T[n];

   // allocate memory
   ia = new long unsigned int [n+1];
   ja = new long unsigned int [nnz];
   a  = new double [nnz];

   Q_lower.makeCompressed();

   for (i = 0; i < n+1; ++i){
      ia[i] = Q_lower.outerIndexPtr()[i]; 
   }  

   for (i = 0; i < nnz; ++i){
      ja[i] = Q_lower.innerIndexPtr()[i];
   }  

   for (i = 0; i < nnz; ++i){
      a[i] = Q_lower.valuePtr()[i];
   }

   // assign b to correct format
   for (int i = 0; i < n*nrhs; i++){
       b[i] = B.data()[i];
       //printf("%f\n", b[i]);
   }

   std::cout << "================== RGF main call ===================" << std::endl;

   RGF<T> *solver;  
   solver = new RGF<T>(ia, ja, a, ns, nt, nb);

   t0 = get_time(0.0);
   solver->factorize();
   t0 = get_time(t0);

   t1 = get_time(0.0);
   solver->solve(x, b, nrhs);
   t1 = get_time(t1);

   t2 = get_time(0.0);
   double logdet = solver->logDet();
   t2 = get_time(t2);

   /*t3 = get_time(0.0);
   solver->RGFdiag(invDiag);
   t3 = get_time(t3);*/

   printf("nnz    : %ld\n", nnz); 
   printf("logdet : %f\n", logdet);
   printf("\n");

   printf("factorize time : %lg\n",t0);
   printf("solve time     : %lg\n",t1);
   printf("total time     : %lg\n", t0+t1);
   //printf("logdet time: %lg\n",t2);
   //printf("partial inverse time: %lg\n",t3);
   printf("\n");

   printf("Residual norm            : %e\n", solver->residualNorm(x, b));
   printf("Residual norm normalized : %e\n", solver->residualNormNormalized(x, b));
   printf("\n");

   /*for (int i = 0; i < nrhs*size; i++)
   {
    printf("%32.24e\n", x[i]);
    //printf("%32.24e\n", invDiag[i]);
   }*/

   // get solution x into matrix X
   MatrixXd X(n, nrhs);
   for (int i = 0; i < n*nrhs; i++){
     X.data()[i] = x[i];
     //printf("%f\n", x[i]);
   }   

   /*MatrixXd res1 = Q*X;
   std::cout << "Q*X:\n" << res1.transpose() << std::endl;*/


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
