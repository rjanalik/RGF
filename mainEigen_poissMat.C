#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#include <armadillo>

#include "read_write_functions.cpp"
#include "generate_matrices/generate_poisson_mat.cpp"
#include "RGF.H"

//#include <likwid.h>


using Eigen::MatrixXd;

typedef Eigen::VectorXd Vect;
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

#define PRINT_MSG

#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif

/* ===================================================================== */

int main(int argc, char* argv[])
{

    if(argc != 1 + 3){
        std::cout << "wrong number of input parameters. " << std::endl;

        std::cerr << "INLA Call : ns nt nb noRhs" << std::endl;

        std::cerr << "[integer:ns]                number of spatial grid points " << std::endl;
        //std::cerr << "[integer:nt]                number of temporal grid points " << std::endl;
        std::cerr << "[integer:nb]                number of fixed effects" << std::endl;
        std::cerr << "[integer:noRhs]             number of right hand sides" << std::endl;

            
        exit(1);
    }

    std::cout << "========= New RGF main call ===========" << std::endl;

    size_t i; // iteration variable

    //std::cout << "reading in example. " << std::endl;

    size_t ns = atoi(argv[1]);
    size_t nt = ns;
    //size_t nt = atoi(argv[2]);
    size_t nb = atoi(argv[2]);
    size_t noRhs = atoi(argv[3]);

    std::cout << "ns = " << ns << ", nt = " << nt << ", nb = " << nb << ", noRhs = " << noRhs << std::endl;

    size_t n  = ns*nt + nb;


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
    std::cout << "Q:\n" << MatrixXd(Q) << std::endl;

    MatrixXi B(n, noRhs);
    generate_int_rhs(B, n, noRhs);

	// =========================================================================== //
	//std::cout << "Converting Eigen Matrices to CSR format. " << std::endl;

	// only take lower triangular part of A
    SpMat Q_lower = Q.triangularView<Lower>(); 
    size_t nnz = Q_lower.nonZeros();
    std::cout << "nnz : " << nnz << std::endl;

    size_t* ia; 
    size_t* ja;
    T* a; 
    T *b;
  	T *x;

  	b  = new T[n*noRhs];
  	x  = new T[n*noRhs];

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

    #if 1
    double t_factorise;
	double t_solve;
	RGF<T> *solver;

	solver = new RGF<T>(ia, ja, a, ns, nt, nb);

	t_factorise = get_time(0.0);
	//solver->solve_equation(GR);
	double flops_factorize = solver->factorize();
	t_factorise = get_time(t_factorise);

	double log_det = solver->logDet();
	printf("logdet: %f\n", log_det);

  	// assign b to correct format
  	for (int i = 0; i < n*noRhs; i++){
	    b[i] = B.data()[i];
	    //printf("%f\n", b[i]);
  	}

  	t_solve = get_time(0.0); 
  	double flops_solve = solver->solve(x, b, noRhs);
  	t_solve = get_time(t_solve);
  	//printf("flops solve:     %f\n", flops_solve);

  	printf("RGF factorise time   : %lg sec\n",t_factorise);
  	printf("RGF solve time       : %lg sec\n",t_solve);

	printf("Residual norm        : %e\n", solver->residualNorm(x, b));
	printf("Res. norm normalised : %e\n", solver->residualNormNormalized(x, b));

    // get solution x into matrix X
    MatrixXd X(n, noRhs);
    for (int i = 0; i < n*noRhs; i++){
        X.data()[i] = x[i];
        //printf("%f\n", b[i]);
    }

    std::cout << "X:\n" << X << std::endl;
    MatrixXd res1 = MatrixXd(Q)*X;
    std::cout << "Q*X:\n" << res1 << std::endl;
    std::cout << "B:\n" << B << std::endl;

    //MatrixXd Bd = MatrixXd(B);
    //MatrixXd res2 = res1 - Bd;
    //std::cout << "Q*X - B = " << (res - B).eval() << std::endl;

  #endif
  
  // free memory
  delete solver;
  delete[] ia;
  delete[] ja;
  delete[] a;
  delete[] b;
  delete[] x;

    
  return 0;


  }
