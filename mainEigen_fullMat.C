#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#include <armadillo>

#include "read_write_functions.cpp"
#include "RGF.H"

//#include <likwid.h>


using Eigen::VectorXd;
using Eigen::MatrixXd;

typedef Eigen::VectorXd Vector;

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

    if(argc != 1 + 4){
        std::cout << "wrong number of input parameters. " << std::endl;

        std::cerr << "INLA Call : ns nt nb path/to/files" << std::endl;

        std::cerr << "[integer:ns]                number of spatial grid points " << std::endl;
        std::cerr << "[integer:nt]                number of temporal grid points " << std::endl;
        std::cerr << "[integer:nb]                number of fixed effects" << std::endl;

        std::cerr << "[string:base_path]          path to folder containing matrix files " << std::endl;
    
        exit(1);
    }

    std::cout << "========= New RGF main call ===========" << std::endl;

    size_t i; // iteration variable

    //std::cout << "reading in example. " << std::endl;

    size_t ns = atoi(argv[1]);
    size_t nt = atoi(argv[2]);
    size_t nb = atoi(argv[3]);
    std::cout << "ns = " << ns << ", nt = " << nt << ", nb = " << nb << std::endl;

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

    std::string base_path = argv[4];    

    /* ---------------- read in matrix ---------------- */

    std::cout << "loading data." << std::endl;
    SpMat Q(n,n);
    std::string Q_file       =  base_path + "/ns" + ns_s + "_nt" + nt_s + "_nb" + nb_s + ".mat";
    std::cout << Q_file << std::endl;
    file_exists(Q_file);
    Q = read_sym_CSR(Q_file);
    //std::cout << "Q = \n" << Q << std::endl;

    Vector rhs(n);
    std::string rhs_file        =  base_path + "/rhs" + to_string(n) + ".txt";
    file_exists(rhs_file);
    rhs = read_matrix(rhs_file, n, 1);

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

  	b  = new T[n];
  	x  = new T[n];

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
  	for (int i = 0; i < n; i++){
	    b[i] = rhs[i];
	    //printf("%f\n", b[i]);
  	}

  	t_solve = get_time(0.0); 
  	double flops_solve = solver->solve(x, b, 1);
  	t_solve = get_time(t_solve);
  	//printf("flops solve:     %f\n", flops_solve);

  	printf("RGF factorise time   : %lg sec\n",t_factorise);
  	printf("RGF solve time       : %lg sec\n",t_solve);

	printf("Residual norm        : %e\n", solver->residualNorm(x, b));
	printf("Res. norm normalised : %e\n", solver->residualNormNormalized(x, b));

  	// create file with solution vector
  	/*std::string sol_x_file_name = "x_sol_RGF_ns" + ns_s + "_nt" + nt_s + "_nb" + nb_s + "_no" + no_s +".dat";
  	std::ofstream sol_x_file(sol_x_file_name,    std::ios::out | std::ios::trunc);

	for (i = 0; i < n; i++) {
		sol_x_file << x[i] << std::endl;
		// sol_x_file << x[i] << std::endl; 
	}

  sol_x_file.close();*/

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
