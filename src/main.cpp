#include <fstream>
/**
 * @file      main.cpp
 * @brief     Runner file that assembles the SPDE and runs the RGF solver
 * @date      Mon Jun 21 08:57:06 2021
 * @author    Radim & Lisa
 *
 */

// main_const_ind

#include "RGF.hpp"
#include "ModelGenerator.hpp"
#include "utilities.h"
#include <iostream>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <armadillo>
#include <omp.h>

using namespace std;

#if 0
typedef CPX T;
    #define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
    #define assign_T(val) val
#endif

void parse_args(int argc, char *argv[]);

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
    size_t ns, nt, nb, no, nu;
    std::string ns_s, nt_s, nb_s, no_s, nu_s;
    arma::vec theta;
    ModelGenerator *model;
    parse_args(argc, argv, base_path, ns, nt, nb, no, theta);
#ifdef DEBUG
    theta.print();
#endif
///////////////////////////////////////////////////////////////////////////////
//                             Generate Model                                //
///////////////////////////////////////////////////////////////////////////////
    model = new ModelGenerator(ns, nt, nb, no, theta, base_path);
    model->construct_model();
    model->assemble_triplet_format();
    printf("\nAll matrices assembled. Passing to RGF solver now.\n");

///////////////////////////////////////////////////////////////////////////////
//                                 RGF SOLVER                                //
///////////////////////////////////////////////////////////////////////////////
    int i;
    double data;
    double t_factorise;
    double t_solve;
    double t_inv;
    size_t n = model->get_n();
    size_t nnz = model->get_nnz();
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

    solver = new RGF<T>(model->triplets.row_idx, model->triplets.col_ptr, model->triplets.val, ns, nt, nb);

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
    printf("Residual norm normalized: %e\n", solver->residualNormNormalized(x, b));

    // create file with solution vector
    std::string sol_x_file_name = base_path + "/x_sol_RGF" + "_ns" + ns_s + "_nt" + nt_s + "_nb" + nb_s + "_no" + no_s + ".dat";
    std::ofstream sol_x_file(sol_x_file_name, std::ios::out | std::ios::trunc);

    for (i = 0; i < n; i++) {
        sol_x_file << x[i] << std::endl;
        // sol_x_file << x[i] << std::endl;
    }
    sol_x_file.close();

    std::string log_file_name = base_path + "/log_RGF_ns" + ns_s + "_nt" + nt_s + "_nb" + nb_s + "_no" + no_s + ".dat";
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
    string sel_inv_file_name =
        base_path + "/RGF_sel_inv_ns" + to_string(ns) + "_nt" + to_string(nt) + "_nb" + nb_s + "_no" + no_s + ".dat";
    std::cout << sel_inv_file_name << endl;
    ofstream sel_inv_file(sel_inv_file_name, ios::out | ::ios::trunc);

    for (int i = 0; i < n; i++) {
        sel_inv_file << invDiag[i] << endl;
    }

    sel_inv_file.close();
    std::cout << "after writing file " << endl;

    // free memory
    delete[] invDiag;
    delete solver;
    delete[] b;
    delete[] x;

#if 0
#endif

    return 0;
}
