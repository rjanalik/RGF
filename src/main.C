#include <math.h>
#include <fstream>
#include <iostream>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include <omp.h>
#include "Utilities.H"

using namespace std;

enum RGF_VERSIONS { BASELINE, ASYNCHRONOUS, ASYNCHRONOUS_2S , BANDED_COPYING, PARDISO_VERSION};
std::map<RGF_VERSIONS, std::string> enum_to_string{{BASELINE, "BASELINE"},
                                                  {ASYNCHRONOUS, "ASYNC"},
                                                  {ASYNCHRONOUS_2S, "ASYNC_TWO_STREAMS"},
                                                  {BANDED_COPYING, "COPYING_NZ_BANDS"},
                                                  {PARDISO_VERSION, "PARDISO"}};
#ifdef BASE
#include "RGF.H"
RGF_VERSIONS rgf_ver = BASELINE;
#elif defined ASYNC
#include "RGF_async.H"
RGF_VERSIONS rgf_ver = ASYNCHRONOUS;
#elif defined ASYNC2S
#include "RGF_async_2s.H"
RGF_VERSIONS rgf_ver = ASYNCHRONOUS_2S;
#elif defined BANDED
#include "RGF_no_zero_copying.H"
RGF_VERSIONS rgf_ver = BANDED_COPYING;
#elif defined PARDISO
#include "PARDISO.H"
RGF_VERSIONS rgf_ver = PARDISO_VERSION;
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);
// #define size_t int
#endif


#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif

int main(int argc, char *argv[])
{
    double overall_time = -omp_get_wtime();
    std::string base_path;
    size_t ns, nt, nb, rows, nnz;
    bool overwrite_results;
    // std::string ns_s, nt_s, nb_s, no_s, nu_s;
    size_t nrhs = 1;
    size_t* ia;
    size_t* ja;
    double* a;
    // arma::vec theta = {};
    // ModelGenerator *model;
    // ALLOCATE MEMORY
    int i;
    double data;
    double t_factorise; double t_solve; double t_inv, t_symbolic_factorise;
#ifndef PARDISO
    // TODO: maybe make nan
    t_symbolic_factorise = 0;
#endif
    T *rhs;
    T *x;
    T *invDiag;
#ifndef PARDISO
    RGF<T> *solver;
#endif
    time_t rawtime;
    parse_args(argc, argv, base_path, ns, nt, nb, overwrite_results);
    rows = ns*nt+nb;
    rhs    = new T[rows];
    x      = new T[rows];
    invDiag= new T[rows];
    std::string mat_path = base_path+"/ns"+std::to_string(ns)+"_nt"+std::to_string(nt)+"_nb"+std::to_string(nb)+".mat";
    std::string rhs_path = base_path+"/rhs"+std::to_string(rows)+".txt";
    utilities::if_not_exists_abort(mat_path);
    utilities::if_not_exists_abort(rhs_path);
    utilities::read_test_matrix_nnz(nnz, mat_path);
    if (rgf_ver == RGF_VERSIONS::BASELINE || rgf_ver == RGF_VERSIONS::PARDISO_VERSION) {
        ia = new size_t[rows + 1];
        ja = new size_t[nnz];
        a = new double[nnz];
    } else {
        cudaMallocHost(&ia, (rows + 1) * sizeof(size_t));
        cudaMallocHost(&ja, nnz * sizeof(size_t));
        cudaMallocHost(&a, nnz * sizeof(size_t));
    }
    // TODO why does read_test_matrix not work?
    // utilities::read_test_matrix(ia, ja, a, rhs, rows, nnz, nrhs, mat_path, rhs_path);
    FILE *F = fopen(mat_path.c_str(), "r");
    double val;
    // read in matrix A, sparse matrix in CSR format
    fscanf(F,"%zu",&rows);
    fscanf(F,"%zu",&rows);
    fscanf(F,"%zu",&nnz);
    for (i = 0; i <= rows; i++){
       fscanf(F,"%zu",&ia[i]);
    }
    for (i = 0; i < ia[rows]; i++){
       fscanf(F,"%zu",&ja[i]);
    }
    for (i = 0; i < ia[rows]; i++){
       fscanf(F,"%lf",&val);
       a[i] = assign_T(val);
    }
    fclose(F);
    // READ RHS
    nrhs   = 1;
    rhs    = new T[nrhs*rows];
    F = fopen(rhs_path.c_str(), "r");
    for (int i = 0; i < nrhs*rows; i++)
    {
        fscanf(F,"%lf",&rhs[i]);
    }
    fclose(F);
#ifdef DEBUG_L2
    utilities::print_header("Inital Matrix Stucture");
    utilities::print_csr(ia, ja, a, rows, nb, true);
    utilities::print_csr(ia, ja, a, rows, nb, false);
#endif
double flops_factorize;
double flops_inv;
double flops_solve;
#ifndef PARDISO
    // load matrix from file
    solver = new RGF<T>(ia, ja, a, ns, nt, nb);
    // FACTORIZATION ===========================================
    t_factorise = get_time(0.0); //solver->solve_equation(GR);
    flops_factorize = solver->factorize();
    t_factorise = get_time(t_factorise);
    double log_det = solver->logDet();
    printf("logdet: %f\n", log_det);
    printf("flops factorize: %f\n", flops_factorize);
    // SOLVE ===========================================
    t_solve = get_time(0.0);
    flops_solve = solver->solve(x, rhs, 1);
    t_solve = get_time(t_solve);
    printf("flops solve:     %f\n", flops_solve);
    // INVERSION ===========================================
    t_inv = get_time(0.0);
    flops_inv = solver->RGFdiag(invDiag);
    t_inv = get_time(t_inv);
    printf("flops inv:      %f\n", flops_inv);
#endif

    // SAVE THE RESULTS ===========================================
    std::ofstream output_fh;
    std::string results_f = "/home/x_pollakgr/RGF/results/tests.csv";
    utilities::if_not_exists_abort(results_f);
    output_fh.open(results_f, overwrite_results ? std::ofstream::trunc : std::ofstream::app);
    output_fh.precision(17);
    std::string sep = "\t";
    if(overwrite_results){
        // WRITE HEADER
        output_fh << "Date/Time"
         << sep << "ns" <<
            sep << "nt" <<
            sep << "nb" <<
            sep << "t_factorize" <<
            sep << "flops_factorize" <<
            sep << "t_solve" <<
            sep << "flops_solve" <<
            sep << "t_inv" <<
            sep << "flops_inv" <<
            sep << "t_symbolic_factorise" <<
            sep << "RGF_Version" << "\n";
    }
    output_fh << std::fixed << utilities::currentDateTime() <<
      sep << std::to_string(ns) <<
      sep << std::to_string(nt) <<
      sep << std::to_string(nb) <<
      sep << t_factorise <<
      sep << flops_factorize <<
      sep << t_solve <<
      sep << flops_solve <<
      sep << t_inv <<
      sep << flops_inv <<
      sep << t_symbolic_factorise <<
      sep << enum_to_string[rgf_ver] <<
      "\n";

    output_fh.close();

#ifdef PARDISO
    printf("PARDISO factorise time: %lg\n",t_factorise);
    printf("PARDISO solve     time: %lg\n",t_solve);
    printf("PARDISO sel inv   time: %lg\n",t_inv);
    printf("PARDISO symbolic factorse time: %lg\n",t_symbolic_factorise);
#else
    printf("RGF factorise time: %lg\n",t_factorise);
    printf("RGF solve     time: %lg\n",t_solve);
    printf("RGF sel inv   time: %lg\n",t_inv);
#endif


#ifndef PARDISO
    printf("Residual norm: %e\n", solver->residualNorm(x, rhs));
    printf("Residual norm normalized: %e\n", solver->residualNormNormalized(x, rhs));
#endif

    // free memory
    if (rgf_ver == RGF_VERSIONS::BASELINE || rgf_ver == RGF_VERSIONS::PARDISO_VERSION) {
        delete[] ia;
        delete[] ja;
        delete[] a;
    } else {
        cudaFreeHost(ia);
        cudaFreeHost(ja);
        cudaFreeHost(a);
    }

    delete[] invDiag;
    delete[] rhs;
    delete[] x;
#ifndef PARDISO
    delete solver;
#endif

    return 0;
}

/************************************************************************************************/
