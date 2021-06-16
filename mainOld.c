#include "RGF.H"
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

/*
Start simulations with RGF NBlock Bmin.dat Bmax.dat M.dat
NBlock: number of blocks of the matrix
Bmin.dat: name of file that contains the first index of each block starting with
0 Bmax.dat: name of file that contains the last index of each block with
Bmax[i]=Bmin[i+1] M.dat: name of file that contains the matrix to work on. Data
stored: size (matrix size) n_nonzeros (number of non-zero elements) fortran
index (0 or 1) index_i index_j real imag (4 columns per matrix entry)
*/

void readCSR(std::string filename, int &n, int &nnz, int *ia, int *ja,
             double *a) {

    fstream fin(filename, ios::in);
    fin >> n;
    fin >> n;
    fin >> nnz;

    // allocate memory
    ia = new int[n + 1];
    ja = new int[nnz];
    a = new double[nnz];

    for (int i = 0; i <= n; i++) {
        fin >> ia[i];
    }

    for (int i = 0; i < ia[n]; i++) {
        fin >> ja[i];
    }

    for (int i = 0; i < ia[n]; i++) {
        fin >> a[i];
    }

    fin.close();
}

arma::sp_mat readCSC(std::string filename) {
    int n_rows;
    int n_cols;
    int nnz;

    ifstream fin(filename);

    fin >> n_rows;
    fin >> n_cols;
    fin >> nnz;

    // allocate memory
    arma::uvec row_ind(nnz);
    arma::uvec col_ptr(n_cols + 1);
    arma::vec a(nnz);

    for (int i = 0; i < nnz; i++) {
        fin >> row_ind[i];
    }

    for (int i = 0; i < n_cols + 1; i++) {
        fin >> col_ptr[i];
        // std::cout <<col_ptr[i] << std::endl;
    }

    for (int i = 0; i < nnz; i++) {
        fin >> a[i];
        // std::cout <<a[i] << std::endl;
    }

    fin.close();

    arma::sp_mat A(row_ind, col_ptr, a, n_rows, n_cols);

    // std::cout << "nonzeros A " << A.n_nonzero << std::endl;
    return A;
}

arma::sp_mat readCSC_sym(std::string filename) {

    int n;
    int nnz;

    fstream fin(filename, ios::in);
    fin >> n;
    fin >> n;
    fin >> nnz;

    // allocate memory
    arma::uvec row_ind(nnz);
    arma::uvec col_ptr(n + 1);
    arma::vec a(nnz);

    for (int i = 0; i < nnz; i++) {
        fin >> row_ind[i];
    }

    for (int i = 0; i < n + 1; i++) {
        fin >> col_ptr[i];
    }

    for (int i = 0; i < nnz; i++) {
        fin >> a[i];
    }

    fin.close();

    arma::sp_mat A_lower(row_ind, col_ptr, a, n, n);
    // create entire matrix
    arma::sp_mat A = arma::symmatl(A_lower);

    return A;
}

arma::mat read_matrix(std::string filename, int n_row, int n_col) {

    arma::mat B(n_row, n_col);
    // std::cout << "size(matrix) : " << arma::size(B) << std::endl;
    B.load(filename,
           arma::raw_ascii); // check later if there is a faster alternative
    // B.print();

    return B;
}

void file_exists(std::string file_name) {
    if (std::fstream{file_name})
        ;
    else {
        std::cerr << file_name
                  << " couldn\'t be opened (not existing or failed to open)\n";
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    cout << "in main" << endl;

    if (argc != 6 + 1) {
        std::cerr << "Pardiso Call : path_to_folder ns nt nb path_to_data no"
                  << std::endl;

        std::cerr << "[string:base_path]          path to folder containing "
                     "spat temp files "
                  << std::endl;
        std::cerr
            << "[integer:ns]                number of spatial grid points "
            << std::endl;
        std::cerr
            << "[integer:nt]                number of temporal grid points "
            << std::endl;
        std::cerr << "[integer:nb]                number of fixed effects"
                  << std::endl;

        std::cerr
            << "[string:base_path_data]          path to folder containing data"
            << std::endl;
        std::cerr << "[integer:no]                number of data samples"
                  << std::endl;

        exit(1);
    }

    // start timer
    double overall_time = -omp_get_wtime();

    std::string base_path = argv[1];
    size_t ns = atoi(argv[2]);
    size_t nt = atoi(argv[3]);
    size_t nb = atoi(argv[4]);

    std::string base_path_data = argv[5];
    size_t no = atoi(argv[6]);

    size_t nu = ns * nt;

    // also save as string
    std::string ns_s = std::to_string(ns);
    std::string nt_s = std::to_string(nt);
    std::string nb_s = std::to_string(nb);
    std::string no_s = std::to_string(no);

    std::string nu_s = std::to_string(nu);

    // check which year's data we are using
    std::string year1 = "2019";
    std::string year2 = "2020";
    std::string year;

    std::size_t found1 = base_path_data.find(year1);
    if (found1 != std::string::npos) {
        std::cout << "found year 2019 in base_path_data " << found1 << '\n';
        year = "2019";
    }

    std::size_t found2 = base_path_data.find(year2);
    if (found2 != std::string::npos) {
        std::cout << "found year 2020 in base_path_data " << found2 << '\n';
        year = "2020";
    }

    std::cout << "./pardiso_call_INLA_sol0_tmp_data " << ns_s << " " << nt_s
              << " " << nb_s << " " << no_s << "\n"
              << std::endl;

    // ------------------- construct file names and check if files exist
    // --------------------------

    // files to construct Q.u depending on HYPERPARAMETERS theta
    std::string const_file =
        base_path + "/const_ns" + ns_s + "_nt" + nt_s + ".dat";
    file_exists(const_file);

    std::string q1s_file = base_path + "/q1s_" + ns_s + ".dat";
    file_exists(q1s_file);
    std::string q2s_file = base_path + "/q2s_" + ns_s + ".dat";
    file_exists(q2s_file);
    std::string q3s_file = base_path + "/q3s_" + ns_s + ".dat";
    file_exists(q3s_file);

    std::string M0_file = base_path + "/M0_" + nt_s + ".dat";
    file_exists(M0_file);
    std::string M1_file = base_path + "/M1_" + nt_s + ".dat";
    file_exists(M1_file);
    std::string M2_file = base_path + "/M2_" + nt_s + ".dat";
    file_exists(M2_file);

    // INDEPENDENT OF HYPERPARAMETERS
    // files containing Q.b
    std::string Qb_file = base_path_data + "/Qb_" + nb_s + "_" + nb_s + ".dat";
    file_exists(Qb_file);

    // files containing B
    std::string B_file = base_path_data + "/B_" + no_s + "_" + nb_s + "_ns" +
                         ns_s + "_nt" + nt_s + ".dat";
    file_exists(B_file);

    // files related to data A.st, y
    std::string A_st_file = base_path_data + "/A_st_" + no_s + "_" +
                            std::to_string(ns * nt) + "_ns" + ns_s + "_nt" +
                            nt_s + ".dat";
    file_exists(A_st_file);
    std::string y_file =
        base_path_data + "/y_" + no_s + "_1_ns" + ns_s + "_nt" + nt_s + ".dat";
    file_exists(y_file);

    // ------------------------------- read in files
    // --------------------------------- //

    // Initialize the random generator
    // arma::arma_rng::set_seed_random();
    arma::vec g(4);

    fstream fin(const_file, ios::in);
    fin >> g[0];
    fin >> g[1];
    fin >> g[2];
    fin >> g[3];
    fin.close();

    // READ IN MATRICES
    // spatial matrices
    arma::sp_mat q1s = readCSC_sym(q1s_file);
    // arma::mat(q1s).print();

    arma::sp_mat q2s = readCSC_sym(q2s_file);
    arma::sp_mat q3s = readCSC_sym(q3s_file);

    // temporal matrices
    arma::sp_mat M0 = readCSC_sym(M0_file);
    arma::sp_mat M1 = readCSC_sym(M1_file);
    arma::sp_mat M2 = readCSC_sym(M2_file);

    // Q.b (sparse & symmetric)
    arma::sp_mat Q_b = readCSC_sym(Qb_file);

    // Qub0  --> create here, all zeros
    arma::sp_mat Qub0(nb, ns * nt);
    Qub0.zeros();

    // Q.e (sparse & symmetric)
    arma::sp_mat Q_e(no, no);
    Q_e = exp(g[0]) * Q_e.eye();

    // A.st (sparse)
    arma::sp_mat A_st = readCSC(A_st_file);
    // std::cout << "non zeros A_st " << A_st.n_nonzero << std::endl;

    // for now B just intercept, all ones
    arma::mat B = read_matrix(B_file, no, nb);

    // y (vector)
    arma::vec y = read_matrix(y_file, no, 1);

    /* ------------- ASSEMBLE MATRICES  ----------------- */

    // Q.u
    arma::sp_mat Q_u(nu, nu);
    Q_u = g[1] * (kron(M0, q3s) + kron(M1 * 2 * g[2], q2s) +
                  kron(M2 * g[2] * g[2], q1s));

    arma::sp_mat Q_u_lower = arma::trimatl(Q_u);
    unsigned int nnz_u = Q_u_lower.n_nonzero;

    // assemble Q.x = block_diagonal(Q.u, Q.b)
    size_t n = size(Q_u)[1] + size(Q_b)[1];
    // std::cout << "n : " << n << std::endl;
    arma::sp_mat Q_x(n, n);
    Q_x(0, 0, size(Q_u)) = Q_u;
    Q_x(0, nu, size(Qub0.t())) = Qub0.t();
    Q_x(nu, 0, size(Qub0)) = Qub0;
    Q_x(nu, nu, size(Q_b)) = Q_b;

    // A.x = cbind(A.st, B)
    arma::sp_mat A_x(no, n);
    A_x(0, 0, size(A_st)) = A_st;
    A_x(0, nu, size(B)) = B;

    // Q.x|y = Q.x + t(A.x), Q.e*A.x
    arma::sp_mat Q_xy(size(Q_x));
    Q_xy = Q_x + A_x.t() * Q_e * A_x;

    arma::vec B_xey = A_x.t() * Q_e * y;

    std::cout << "size(Q_xy)" << size(Q_xy) << std::endl;
    std::cout << "size(B_xey)" << size(B_xey) << std::endl;

    // TAKE ENTIRE MATRIX FOR THIS SOLVER
    arma::sp_mat Q_xy_lower = arma::trimatl(Q_xy);

    // require CSR format
    size_t nnz = Q_xy_lower.n_nonzero;

    std::cout << "number of non zeros : " << nnz << std::endl;

    size_t *ia;
    size_t *ja;
    double *a;

    // allocate memory
    ia = new size_t[n + 1];
    ja = new size_t[nnz];
    a = new double[nnz];

    for (int i = 0; i < n + 1; ++i) {
        ia[i] = Q_xy_lower.col_ptrs[i];
    }

    for (int i = 0; i < nnz; ++i) {
        ja[i] = Q_xy_lower.row_indices[i];
    }

    for (int i = 0; i < nnz; ++i) {
        a[i] = Q_xy_lower.values[i];
    }

    printf("\nAll matrices assembled. Passing to RGF solver now.\n");

    // -------------------------------------------------------------------------------------------
    // //
    // -------------------------------------- RGF SOLVER
    // ----------------------------------------- //
    // -------------------------------------------------------------------------------------------
    // //

    // FILE *F1,*F2;
    int i;
    double data;
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

    // assign b to correct format
    for (int i = 0; i < n; i++) {
        b[i] = B_xey[i];
        // printf("%f\n", b[i]);
    }

    solver = new RGF<T>(ia, ja, a, ns, nt, nb);

    t_factorise = get_time(0.0);
    // solver->solve_equation(GR);
    double flops_factorize = solver->factorize();
    t_factorise = get_time(t_factorise);

    double log_det = solver->logDet();
    printf("logdet: %f\n", log_det);

    printf("flops factorize: %f\n", flops_factorize);

    // write this to file
    /*std::string L_factor_file_name = base_path + "/L_factor_RGF"  + "_ns" +
    ns_s + "_nt" + nt_s + "_nb" + nb_s + "_no" + no_s + ".dat"; std::ofstream
    L_factor_file(L_factor_file_name,    std::ios::out | std::ios::trunc);

    L_factor_file << n << std::endl;
    L_factor_file << n << std::endl;
    L_factor_file << M->n_nonzeros << std::endl;

    for (int i = 0; i < M->size+1; ++i){
      L_factor_file << M->edge_i[i] << std::endl;
    }
     for (int i = 0; i < M->n_nonzeros; ++i){
      L_factor_file << M->index_j[i] << std::endl;
    }
    for (int i = 0; i < M->n_nonzeros; ++i){
      L_factor_file << MF[i] << std::endl;
    }

    L_factor_file.close(); */

    t_solve = get_time(0.0);
    double flops_solve = solver->solve(x, b, 1);
    t_solve = get_time(t_solve);
    printf("flops solve:     %f\n", flops_solve);

    t_inv = get_time(0.0);
    double flops_inv = solver->RGFdiag(invDiag);
    t_inv = get_time(t_inv);
    printf("flops inv:      %f\n", flops_inv);

    // exit(1);

    printf("RGF factorise time: %lg\n", t_factorise);
    printf("RGF solve     time: %lg\n", t_solve);
    printf("RGF sel inv   time: %lg\n", t_inv);

    printf("Residual norm: %e\n", solver->residualNorm(x, b));
    printf("Residual norm normalized: %e\n",
           solver->residualNormNormalized(x, b));

    /*for (int i = 0; i < nrhs*ns*nt; i++){
      printf("x[%d] = %f\n", i, b[i]);
    }*/

    // create file with solution vector
    std::string sol_x_file_name = base_path_data + "/x_sol_RGF" + "_ns" + ns_s +
                                  "_nt" + nt_s + "_nb" + nb_s + "_no" + no_s +
                                  "_year" + year + ".dat";
    std::ofstream sol_x_file(sol_x_file_name, std::ios::out | std::ios::trunc);

    for (i = 0; i < n; i++) {
        sol_x_file << x[i] << std::endl;
        // sol_x_file << x[i] << std::endl;
    }

    sol_x_file.close();

    std::string log_file_name = base_path_data + "/log_RGF_ns" + ns_s + "_nt" +
                                nt_s + "_nb" + nb_s + "_no" + no_s + "_year" +
                                year + ".dat";
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
    string sel_inv_file_name = base_path_data + "/RGF_sel_inv_ns" +
                               to_string(ns) + "_nt" + to_string(nt) + "_nb" +
                               nb_s + "_no" + no_s + ".dat";
    cout << sel_inv_file_name << endl;
    ofstream sel_inv_file(sel_inv_file_name, ios::out | ::ios::trunc);

    for (int i = 0; i < n; i++) {
        sel_inv_file << invDiag[i] << endl;
    }

    sel_inv_file.close();
    cout << "after writing file " << endl;

    // free memory
    delete[] invDiag;
    delete solver;
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
    delete[] x;

    return 0;
}

/************************************************************************************************/
