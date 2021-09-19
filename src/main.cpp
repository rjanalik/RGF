#include <fstream>
/**
 * @file      main.cpp
 * @brief     Runner file that assembles the SPDE and runs the RGF solver
 * @date      Mon Jun 21 08:57:06 2021
 * @author    Radim & Lisa
 *
 */

// main_const_ind

#include "ModelGenerator.hpp"
#include "RGF.hpp"
#include "utilities.hpp"
#include <iostream>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <armadillo>
#include <omp.h>

#ifdef LISA_VERSION
arma::sp_mat construct_Q_spatial(arma::vec theta_u, \
                         arma::sp_mat c0, arma::sp_mat g1, arma::sp_mat g2){

  arma::sp_mat Qs(c0.n_rows, c0.n_cols);
  // Qs <- g[1]^2*Qgk.fun(sfem, g[2], order)
  // return(g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2)
  Qs = pow(theta_u[0],2)*(pow(theta_u[1], 4) * c0 + 2*pow(theta_u[1],2) * g1 + g2);
  // extract triplet indices and insert into Qx

  return Qs;

}

// SPDE discretisation -- matrix construction
//void construct_Q_spatial(arma::vec& theta, arma::sp_mat* Qs, \
//                         arma::sp_mat* c0, arma::sp_mat* g1, arma::sp_mat* g2){
arma::sp_mat construct_Q_spat_temp(arma::vec theta_u, \
                         arma::sp_mat c0, arma::sp_mat g1, arma::sp_mat g2, arma::sp_mat g3, \
                         arma::sp_mat M0, arma::sp_mat M1, arma::sp_mat M2){

  int n_st = c0.n_rows * M0.n_rows;
  arma::sp_mat Qst(n_st, n_st);

  // g^2 * fem$c0 + fem$g1
  arma::sp_mat q1s = pow(theta_u[2],2) * c0 + g1;
  //std::cout << "q1s :" << std::endl;
  //arma::mat(q1s).submat(0,0,10,10).print();

  // g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2
  arma::sp_mat q2s = pow(theta_u[2], 4) * c0 + 2 * pow(theta_u[2],2) * g1 + g2;

  // g^6 * fem$c0 + 3 * g^4 * fem$g1 + 3 * g^2 * fem$g2 + fem$g3
  arma::sp_mat q3s = pow(theta_u[2], 6) * c0 + 3 * pow(theta_u[2],4) * g1 + 3 * pow(theta_u[2],2) * g2 + g3;

  // assemble overall precision matrix Q.st
  Qst = theta_u[0]*(kron(M0, q3s) + kron(M1 * 2 * theta_u[1], q2s) +  kron(M2 * theta_u[1] * theta_u[1], q1s));

  return Qst;

}
#endif

using namespace std;

#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif

/**
 * @brief Runs Recursive Green Function algorithm
 * @details Solve block diagonal arrow matrices using the RGF algorithm
 * @param[in] argc: Argument counter
 * @param[out] argv: Command line arguments:
 * - [1]: path to data folder - \p base_path
 * - [2]: number of spatial grid points - \p ns
 * - [3]: number of temporal grid points - \p nt
 * - [4]: number of fixed effects - \p nd/nb
 * - [5]: number of data samples - \p no
 * @return saves the selective inverse to a file
 */
int main(int argc, char *argv[]) {

  // start timer
  double overall_time = -omp_get_wtime();
  std::string base_path;
  size_t ns, nt, nb, no;
  std::string ns_s, nt_s, nb_s, no_s, nu_s;
  arma::vec theta = {};
  ModelGenerator *model;
  parse_args(argc, argv, base_path, ns, nt, nb, no);
  ///////////////////////////////////////////////////////////////////////////////
  //                             Generate Model //
  ///////////////////////////////////////////////////////////////////////////////
  printf("\nStart Assembling matrices.\n");
  model = new ModelGenerator(ns, nt, nb, no, theta, base_path);
#ifdef DEBUG
  // model->get_theta().print();
#endif
#ifdef LISA_VERSION
  // spatial matrices
  arma::sp_mat c0 = model->c0_;
  arma::sp_mat g1 = model->g1_;
  arma::sp_mat g2 = model->g2_;

  arma::sp_mat g3 = model->g3_;
  arma::sp_mat M0 = model->M0_;
  arma::sp_mat M1 = model->M1_;
  arma::sp_mat M2 = model->M2_;
  arma::sp_mat Qu;

  arma::sp_mat Ax = model->Ax_;
  arma::vec y = model->y_;
  size_t nu                = ns*nt;
  theta = model->get_theta();

 // assemble Qs for given theta
  if(nt == 1){

    Qu = construct_Q_spatial(exp(theta(arma::span(1,2))), c0, g1, g2);
    // arma::mat(Qu).submat(0,0,10,10).print();

  } else {

    std::cout << exp(theta(arma::span(1,3))) << std::endl;
    Qu = construct_Q_spat_temp(exp(theta(arma::span(1,3))), c0, g1, g2, g3, M0, M1, M2);
    // arma::mat(Qu).submat(0,0,10,10).print();

  }

  // Q.b : diag(1e-5), dim : nb nb
  arma::sp_mat Qb = 1e-5 * arma::speye(nb, nb);
  //Qb.print();

// void ModelGenerator::construct_Qxy_lower() {
  // Qub0  --> create here, all zeros
  arma::sp_mat Qub0(nb, ns*nt);
  Qub0.zeros();

  /* ------------- ASSEMBLE MATRICES  ----------------- */

  // assemble Q.x = block_diagonal(Q.u, Q.b)
  size_t n = size(Qu)[1] + size(Qb)[1];
  std::cout << "n : " << n << std::endl;
  arma::sp_mat Qx(n,n);
  Qx(0,0, size(Qu))          = Qu;
  Qx(0,nu, size(Qub0.t()))    = Qub0.t();
  Qx(nu, 0, size(Qub0))       = Qub0;
  Qx(nu, nu, size(Qb))       = Qb;
  //arma::mat(Qx).submat(0,0,10,10).print();

  // Q.x|y = Q.x + t(A.x), Q.e*A.x
  arma::sp_mat Qxy(size(Qx));
  Qxy = Qx + exp(theta[0])*Ax.t()*Ax;
  //arma::mat(Qx).submat(0,0,10,10).print();


  arma::vec bxy = exp(theta[0])*Ax.t()*y;
  //bxy.print();

  //std::cout << "size(Qxy)" << size(Qxy) << std::endl;
  //std::cout << "size(bxy)" << size(bxy) << std::endl;

  // TAKE ENTIRE MATRIX FOR THIS SOLVER
  arma::sp_mat Qxy_lower = arma::trimatl(Qxy);

  // require CSR format
  size_t nnz = Qxy_lower.n_nonzero;
#else
  model->construct_model();
  model->assemble_triplet_format();
  size_t n = model->get_n();
  size_t nnz = model->get_nnz();
#endif
  ///////////////////////////////////////////////////////////////////////////////
  //                                 RGF SOLVER //
  ///////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  printf("\nAll matrices assembled. Passing to RGF solver now.\n");
#endif
  size_t i;
  double t_factorise;
  double t_solve;
  double t_inv;
  T *b;
  T *x;
  T *invDiag;
  RGF<T> *solver;
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  printf("The current date/time is: %s\n", asctime(timeinfo));

  b = new T[n];
  x = new T[n];
  invDiag = new T[n];

// Initialize Solver ///////////////////////////////////////////////////////////
#ifdef LISA_VERSION
  size_t* ia;
  size_t* ja;
  double* a;

  // allocate memory
  ia = new size_t [n+1];
  ja = new size_t [nnz];
  a = new double [nnz];
  for (int i = 0; i < nnz; ++i){
    ja[i] = Qxy_lower.row_indices[i];
    //std::cout << ja[i] << std::endl;

  }

  for (int i = 0; i < n+1; ++i){
    ia[i] = Qxy_lower.col_ptrs[i];
    //std::cout << ia[i] << std::endl;
  }

  for (int i = 0; i < nnz; ++i){
    a[i] = Qxy_lower.values[i];
    //std::cout << a[i] << std::endl;

  }
  solver = new RGF<T>(ia, ja, a, ns, nt, nb);
#else
  solver = new RGF<T>(model->triplets.ia, model->triplets.ja, model->triplets.a, ns, nt, nb);
#endif

  t_factorise = get_time(0.0);
  // solver->solve_equation(GR);
  double flops_factorize = solver->factorize();
  t_factorise = get_time(t_factorise);

  double log_det = solver->logDet();
  printf("logdet: %f\n", log_det);

  printf("flops factorize: %f\n", flops_factorize);

  t_solve = get_time(0.0);
  double flops_solve = solver->solve(x, b, 1);
  t_solve = get_time(t_solve);
  printf("flops solve:     %f\n", flops_solve);

  t_inv = get_time(0.0);
  double flops_inv = solver->RGFdiag(invDiag);
  t_inv = get_time(t_inv);
  printf("flops inv:      %f\n", flops_inv);

  printf("RGF factorise time: %lg\n", t_factorise);
  printf("RGF solve     time: %lg\n", t_solve);
  printf("RGF sel inv   time: %lg\n", t_inv);

  printf("Residual norm: %e\n", solver->residualNorm(x, b));
  printf("Residual norm normalized: %e\n",
         solver->residualNormNormalized(x, b));

  // create file with solution vector
  std::string sol_x_file_name = base_path + "/x_sol_RGF" + "_ns" + ns_s +
                                "_nt" + nt_s + "_nb" + nb_s + "_no" + no_s +
                                ".dat";
  std::ofstream sol_x_file(sol_x_file_name, std::ios::out | std::ios::trunc);

  for (i = 0; i < n; i++) {
    sol_x_file << x[i] << std::endl;
    // sol_x_file << x[i] << std::endl;
  }
  sol_x_file.close();

  std::string log_file_name = base_path + "/log_RGF_ns" + ns_s + "_nt" + nt_s +
                              "_nb" + nb_s + "_no" + no_s + ".dat";
  std::ofstream log_file(log_file_name);
  log_file << ns << std::endl;
  log_file << nt << std::endl;
  log_file << nb << std::endl;
  log_file << no << std::endl;
  log_file << n << std::endl;
  log_file << nnz << std::endl;
  log_file << "RGF" << std::endl;
  log_file << log_det << std::endl;
  log_file << "0.0" << std::endl;
  log_file << t_factorise << std::endl;
  log_file << t_solve << std::endl;
  log_file << t_inv << std::endl;
  log_file << flops_factorize << std::endl;
  log_file << flops_solve << std::endl;
  log_file << flops_inv << std::endl;

  log_file.close();

  // print/write diag
  string sel_inv_file_name = base_path + "/RGF_sel_inv_ns" + to_string(ns) +
                             "_nt" + to_string(nt) + "_nb" + nb_s + "_no" +
                             no_s + ".dat";
  std::cout << sel_inv_file_name << endl;
  ofstream sel_inv_file(sel_inv_file_name, ios::out | ::ios::trunc);

  for (size_t i = 0; i < n; i++) {
    sel_inv_file << invDiag[i] << endl;
  }

  sel_inv_file.close();
  std::cout << "after writing file " << endl;

#ifdef LISA_VERSION
  delete[] a;
  delete[] ia;
  delete[] ja;
#else
  // free memory
  delete model;
  delete solver;
#endif
  delete[] b;
  delete[] x;
  delete[] invDiag;
  return 0;
#if 0
#endif
}
