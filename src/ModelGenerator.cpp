#include "ModelGenerator.hpp"
#include <string>

using namespace utilities;

ModelGenerator::ModelGenerator(size_t ns, size_t nt, size_t nb, size_t no, arma::vec theta, std::string base_path)
    : ns_(ns), nt_(nt), nb_(nb), no_(no), theta_(theta), base_path_(base_path) {
    nu_ = ns_ * nt_;

    if(theta_.is_empty()){
        if(nt_ == 1){
            theta_ = {-1.5,-5,-2};
        } else {
            theta_ = {5, -10, 2.5, 1};
        }
    }
    getData();
};
///////////////////////////////////////////////////////////////////////////////
//                          Model creation routines                          //
///////////////////////////////////////////////////////////////////////////////
void ModelGenerator::construct_model() {
    print_header();
    construct_Qu();
    construct_Qb();
    construct_Qxy_lower();
    construct_b();
}

/**
 * @brief Constructs Qu
 * @details selects between:
   - Qu_spatial
   - Qu_spatiotemporal
 */
void ModelGenerator::construct_Qu() {
    if (nt_ == 1){
#ifdef DEBUG
        print_header("1. construct_Qu_spatial(exp(theta_(arma::span(1, 2))))");
#endif
        construct_Qu_spatial(exp(theta_(arma::span(1, 2))));
    }
    else{
#ifdef DEBUG
        print_header("1. construct_Qu()construct_Qu_spatio_temporal(exp(theta_(arma::span(1, 3))))");
#endif
        construct_Qu_spatio_temporal(exp(theta_(arma::span(1, 3))));
    }
}

/**
 * @brief Summary
 * @details Description
 * @param[in] theta_u Description
 * @param[out] Qu g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2
 */
void ModelGenerator::construct_Qu_spatial(arma::vec theta_u) {
    // Qs <- g[1]^2*Qgk.fun(sfem, g[2], order)
    Qu.resize(c0.n_rows, c0.n_cols);
    Qu = pow(theta_u[0], 2) * (pow(theta_u[1], 4) * c0 + 2 * pow(theta_u[1], 2) * g1 + g2);
}

void ModelGenerator::construct_Qu_spatio_temporal(arma::vec theta_u) {
    int n_st = c0.n_rows * M0.n_rows;
    Qu.resize(n_st, n_st);
    // g^2 * fem$c0 + fem$g1
    arma::sp_mat q1s = pow(theta_u[2], 2) * c0 + g1;
    // std::cout << "q1s :" << std::endl;
    // arma::mat(q1s).submat(0,0,10,10).print();
    // g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2
    arma::sp_mat q2s = pow(theta_u[2], 4) * c0 + 2 * pow(theta_u[2], 2) * g1 + g2;
    // g^6 * fem$c0 + 3 * g^4 * fem$g1 + 3 * g^2 * fem$g2 + fem$g3
    arma::sp_mat q3s = pow(theta_u[2], 6) * c0 + 3 * pow(theta_u[2], 4) * g1 + 3 * pow(theta_u[2], 2) * g2 + g3;
    // assemble overall precision matrix Q.st
    Qu = theta_u[0] * (kron(M0, q3s) + kron(M1 * 2 * theta_u[1], q2s) + kron(M2 * theta_u[1] * theta_u[1], q1s));
}
void ModelGenerator::construct_Qb() {
#ifdef DEBUG
    print_header("2. construct_Qb()");
#endif
    // Q.b : diag(1e-5), dim : nb nb
    Qb = 1e-5 * arma::speye(nb_, nb_);
#ifdef DEBUG
        Qb.print();
#endif
}

/**
 * @brief Summary
 * @details Description
 */
void ModelGenerator::construct_Qxy_lower() {
#ifdef DEBUG
    print_header("3. ModelGenerator::construct_Qxy_lower()");
#endif
    // Qub0  --> create here, all zeros
    arma::sp_mat Qub0(nb_, ns_ * nt_);
    Qub0.zeros();
    // assemble Q.x = block_diagonal(Q.u, Q.b)
    // n = ns+nb i.e. spatial + fixed effects
    // TODO Check if this is correct with LISA
    n_ = size(Qu)[1] + size(Qb)[1];
    // std::cout << "n : " << n << std::endl;
    arma::sp_mat Qx(n_, n_);
    Qx(0, 0, size(Qu)) = Qu;
#ifdef DEBUG
    std::cout << "n_ = " << n_ << std::endl;
    std::cout << "nu_ = " << nu_ << std::endl;
    std::cout << "size(Qu)[1] = " << size(Qu)[1] << std::endl;
    std::cout << "size(Qu)[2] = " << size(Qu)[2] << std::endl;
    std::cout << "size(Qb)[1] = " << size(Qb)[1] << std::endl;
    std::cout << "size(Qb)[2] = " << size(Qb)[2] << std::endl;
    std::cout << "size(Qub0)) = " << size(Qub0) << std::endl;
    std::cout << "size(Qub0.t())) = " << size(Qub0.t()) << std::endl;
#endif
    Qx(0, nu_, size(Qub0.t())) = Qub0.t();
    Qx(nu_, 0, size(Qub0)) = Qub0;
    Qx(nu_, nu_, size(Qb)) = Qb;
    // arma::mat(Qx).submat(0,0,10,10).print();
    // Q.x|y = Q.x + t(A.x), Q.e*A.x
    arma::sp_mat Qxy(size(Qx));
    Qxy = Qx + exp(theta_[0]) * Ax.t() * Ax;
    // arma::mat(Qx).submat(0,0,10,10).print();
    Qxy_lower = arma::trimatl(Qxy);
    nnz_ = Qxy_lower.n_nonzero;
}

void ModelGenerator::construct_b() {
#ifdef DEBUG
    print_header("4. construct_b()");
#endif
    b = exp(theta_[0]) * Ax.t() * y;
#ifdef DEBUG
    b.print();
#endif
}

void ModelGenerator::assemble_triplet_format(){
    triplets.ia = new size_t[n_+1];
    triplets.ja = new size_t[nnz_];
    triplets.a = new double[nnz_];
    #ifdef DEBUG
    print_header("Assembling Triplets");
    #endif
    for (size_t i = 0; i < nnz_; ++i){
        triplets.ja[i] = Qxy_lower.row_indices[i];
        #ifdef DEBUG
        // std::cout << triplets.ja[i] << std::endl;
        #endif
    }
    for (size_t i = 0; i < n_+1; ++i){
        triplets.ia[i] = Qxy_lower.col_ptrs[i];
        #ifdef DEBUG
        // std::cout << triplets.col_ptr[i] << std::endl;
        #endif
    }
    for (size_t i = 0; i < nnz_; ++i){
        triplets.a[i] = Qxy_lower.values[i];
        #ifdef DEBUG
        // std::cout << triplets.a[i] << std::endl;
        #endif
    }
    #ifdef DEBUG
    std::cout << "n = " << n_ << std::endl;
    std::cout << "nnz = " << nnz_ << std::endl;
    print_header();
    #endif
}
///////////////////////////////////////////////////////////////////////////////
//                              Getter & Setter                              //
///////////////////////////////////////////////////////////////////////////////
size_t ModelGenerator::get_n(){
    return n_;
}
size_t ModelGenerator::get_nnz(){
    return nnz_;
}
///////////////////////////////////////////////////////////////////////////////
//                               Data Fetching                               //
///////////////////////////////////////////////////////////////////////////////
void ModelGenerator::getData() {
    // spatial matrices
    getSpatialData();
    if (nt_ > 1)
        getTemporalData();
};
void ModelGenerator::getSpatialData() {
    std::string ns_s = std::to_string(ns_);
    std::string nt_s = std::to_string(nt_);
    std::string no_s = std::to_string(no_);

    std::string c0_file = base_path_ + "/c0_" + ns_s + ".dat";
    std::string g1_file = base_path_ + "/g1_" + ns_s + ".dat";
    std::string g2_file = base_path_ + "/g2_" + ns_s + ".dat";
    std::string filename_ending = "_ns" + ns_s + "_nt" + nt_s + ".dat";
    std::string Ax_file = base_path_ + "/Ax_" + no_s + "_" + std::to_string(nu_ + nb_) + filename_ending;
    std::string y_file = base_path_ + "/y_" + no_s + "_1" + filename_ending;
    if_not_exists_abort({c0_file, g1_file, g2_file, Ax_file, y_file});
    c0 = readCSC_sym(c0_file);
    g1 = readCSC_sym(g1_file);
    g2 = readCSC_sym(g2_file);
    Ax = readCSC(Ax_file);
    y = read_matrix(y_file, no_, 1);
};
void ModelGenerator::getTemporalData() {
    std::string ns_s = std::to_string(ns_);
    std::string nt_s = std::to_string(nt_);

    std::string g3_file = base_path_ + "/g3_" + ns_s + ".dat";
    std::string M0_file = base_path_ + "/M0_" + nt_s + ".dat";
    std::string M1_file = base_path_ + "/M1_" + nt_s + ".dat";
    std::string M2_file = base_path_ + "/M2_" + nt_s + ".dat";
    if_not_exists_abort({g3_file, M0_file, M1_file, M2_file});
    g3 = readCSC_sym(g3_file);
    M0 = readCSC_sym(M0_file);
    M1 = readCSC_sym(M1_file);
    M2 = readCSC_sym(M2_file);
#ifdef DEBUG
    arma::mat(g3).submat(0, 0, 10, 10).print();
    arma::mat(M0).submat(0, 0, nt_ - 1, nt_ - 1).print();
    arma::mat(M1).submat(0, 0, nt_ - 1, nt_ - 1).print();
    arma::mat(M2).submat(0, 0, nt_ - 1, nt_ - 1).print();
#endif
};

///////////////////////////////////////////////////////////////////////////////
//                             Utility Functions                             //
///////////////////////////////////////////////////////////////////////////////

bool ModelGenerator::file_exists(const std::string &file_name) { return std::fstream{file_name} ? true : false; }
void ModelGenerator::if_not_exists_abort(std::string const file_name) {
    if(file_exists(file_name))
        return;
    std::cerr << file_name
              << " couldn\'t be opened (not existing or failed to open)\n";
    exit(1);
}
void ModelGenerator::if_not_exists_abort(std::initializer_list<std::string> const file_names) {
    for (std::string file_name : file_names)
        if_not_exists_abort(file_name);
}

arma::mat ModelGenerator::read_matrix(std::string filename, int n_row, int n_col) {
    arma::mat B(n_row, n_col);
    // std::cout << "size(matrix) : " << arma::size(B) << std::endl;
    B.load(filename, arma::raw_ascii); // check later if there is a faster alternative
    // B.print();
    return B;
}

arma::sp_mat ModelGenerator::readCSC(std::string filename){
  int n_rows, n_cols, nnz;
  std::ifstream fin(filename);
  fin >> n_rows;
  fin >> n_cols;
  fin >> nnz;
   // allocate memory
  arma::uvec row_ind(nnz);
  arma::uvec col_ptr(n_cols+1);
  arma::vec a(nnz);
  for (int i = 0; i < nnz; i++)
      fin >> row_ind[i];
  for (int i = 0; i < n_cols+1; i++)
      fin >> col_ptr[i];
  for (int i = 0; i < nnz; i++)
      fin >> a[i];
  fin.close();
  arma::sp_mat A(row_ind, col_ptr, a, n_rows, n_cols);
  //std::cout << "nonzeros A " << A.n_nonzero << std::endl;
  return A;
}

arma::sp_mat ModelGenerator::readCSC_sym(std::string filename) {
    std::ifstream fin(filename, std::ios::in);
    size_t n, nnz;
    fin >> n;
    fin >> n;
    fin >> nnz;
    arma::uvec row_ind(nnz);
    arma::uvec col_ptr(n + 1);
    arma::vec a(nnz);

    for (size_t i = 0; i < nnz; i++)
        fin >> row_ind[i];
    for (size_t i = 0; i < n + 1; i++)
        fin >> col_ptr[i];
    for (size_t i = 0; i < nnz; i++)
        fin >> a[i];
    fin.close();
    arma::sp_mat A_lower(row_ind, col_ptr, a, n, n);
    // create entire matrix
    arma::sp_mat A = arma::symmatl(A_lower);

    return A;
}
