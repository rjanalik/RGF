// generate_poission_mat.cpp

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

typedef Eigen::VectorXd Vect;
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

using namespace Eigen;


void generate_adapt_poisson_mat(SpMat& Q, int ns, int nb){

	int n = ns*ns + nb;

	// set up off-diagonals -1's (TODOL: problem with transpose ... therefore upper & lower ... )
	SpMat U(ns,ns);
	U.reserve(ns-1);

	for(int j=0; j<ns-1; ++j)
	    U.insert(j,j+1) = -1;

	SpMat Ut(ns,ns);
	Ut.reserve(ns-1);

	for(int j=0; j<ns-1; ++j)
	    U.insert(j+1,j) = -1;

	// initialize identity 
	SpMat D(ns,ns);
	D.setIdentity();

	// assemble diagonal block
	SpMat A_ns = U + Ut + 2*D;

	//std::cout << "Q:\n" << A_ns << std::endl;

	// assemble using kronecker products (requires ns=nt)
	//A = kron(A_ns, I) + kron(I, A_ns);
	SpMat A = KroneckerProductSparse<SpMat, SpMat>(A_ns, D) + KroneckerProductSparse<SpMat, SpMat>(D, A_ns);
	A = A + nb*MatrixXd::Identity(ns*ns, ns*ns); // to ensure positive definiteness !!

	// add random "fixed effects"
	MatrixXd randMat = 0.5*(MatrixXd::Random(nb, ns*ns)-MatrixXd::Ones(nb, ns*ns));
	//std::cout << "randMat:\n" << randMat << std::endl;

	MatrixXd randDiagBl =  0.5*(MatrixXd::Random(nb, nb)-MatrixXd::Ones(nb, nb)) + (n+1)*MatrixXd::Identity(nb,nb);
	randDiagBl.triangularView<Upper>() = randDiagBl.transpose(); // make this block symmetric!!
	//std::cout << "randDiagBl:\n" << randDiagBl << std::endl;

	// assemble everything together .. make everything dense to make it easier ....
	MatrixXd Q_dense(n,n);
	Q_dense << MatrixXd(A), randMat.transpose(), randMat, randDiagBl;
	//std::cout << "Q dense :\n" << Q_dense << std::endl; 

	//std::cout << "eigenvalues :" << Q_dense.eigenvalues().transpose() << std::endl;

	Q = Q_dense.sparseView();

} 


void generate_int_rhs(MatrixXi& B, int n, int nb){

  	MatrixXd temp = (MatrixXd::Random(n,nb)+MatrixXd::Ones(n,nb))*5;  
  	B = temp.cast<int>();    
  	//std::cout << "B:\n" << B << std::endl;
}


#if 0
int main(int argc, char* argv[]){

    if(argc != 1 + 2){
        std::cout << "wrong number of input parameters. " << std::endl;

        std::cerr << "[integer:ns]                number of spatial grid points " << std::endl;
        std::cerr << "[integer:nt]                number of temporal grid points " << std::endl;
        std::cerr << "[integer:nb]                number of fixed effects" << std::endl;

        exit(1);
    }

    size_t ns = atoi(argv[1]);
    size_t nt = ns; // at the moment ns = nt ... 
    //size_t nt = atoi(argv[2]);
    size_t nb = atoi(argv[2]);
    std::cout << "ns = " << ns << ", nt = " << nt << ", nb = " << nb << std::endl;
    int n = ns*nt + nb;

    Vect b(n);
    SpMat Q(n,n);    

    generate_adapt_poisson_mat(Q, ns, nb);
    std::cout << "Q =\n" << Q << std::endl;



	return 0;
}

#endif



