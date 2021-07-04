#include "ModelGenerator.hpp"
#include <string>

ModelGenerator::ModelGenerator(size_t ns, size_t nt, size_t nb, size_t no, arma::vec theta, std::string base_path)
    : ns_(ns), nt_(nt), nb_(nb), no_(no), theta_(theta), base_path_(base_path) {
    nu_ = ns_ * nt_;
    ns_s = std::to_string(ns_);
    nt_s = std::to_string(nt_);
    nb_s = std::to_string(nb_);
    no_s = std::to_string(no_);
    nu_s = std::to_string(nu_);
};
///////////////////////////////////////////////////////////////////////////////
//                          Model creation routines                          //
///////////////////////////////////////////////////////////////////////////////
void ModelGenerator::construct_model() {
    construct_Qu();
    construct_Qb();
    construct_Qxy_lower();
    construct_b();
}
void ModelGenerator::construct_Qu() {
    if (nt_ > 1)
        construct_Qu_spatial(exp(theta_(arma::span(1, 2))));
    else
        construct_Qu_spatio_temporal(exp(theta_(arma::span(1, 3))));
}

/**
 * @brief Summary
 * @details Description
 * @param[in] theta_u Description
 * @param[out] Qu g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2
 */
void ModelGenerator::construct_Qu_spatial(arma::vec theta_u) {
    // Qs <- g[1]^2*Qgk.fun(sfem, g[2], order)
    arma::sp_mat Qs(data.c0.n_rows, data.c0.n_cols);
    data.Qu = pow(theta_u[0], 2) * (pow(theta_u[1], 4) * data.c0 + 2 * pow(theta_u[1], 2) * data.g1 + data.g2);
}

void ModelGenerator::construct_Qu_spatio_temporal(arma::vec theta_u) {
    int n_st = data.c0.n_rows * data.M0.n_rows;
    arma::sp_mat Qst(n_st, n_st);
    // g^2 * fem$c0 + fem$g1
    arma::sp_mat q1s = pow(theta_u[2], 2) * data.c0 + data.g1;
    // std::cout << "q1s :" << std::endl;
    // arma::mat(q1s).submat(0,0,10,10).print();
    // g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2
    arma::sp_mat q2s = pow(theta_u[2], 4) * data.c0 + 2 * pow(theta_u[2], 2) * data.g1 + data.g2;
    // g^6 * fem$c0 + 3 * g^4 * fem$g1 + 3 * g^2 * fem$g2 + fem$g3
    arma::sp_mat q3s = pow(theta_u[2], 6) * data.c0 + 3 * pow(theta_u[2], 4) * data.g1 + 3 * pow(theta_u[2], 2) * data.g2 + data.g3;
    // assemble overall precision matrix Q.st
    data.Qu =
        theta_u[0] * (kron(data.M0, q3s) + kron(data.M1 * 2 * theta_u[1], q2s) + kron(data.M2 * theta_u[1] * theta_u[1], q1s));
}
void ModelGenerator::construct_Qb() {
    // Q.b : diag(1e-5), dim : nb nb
    data.Qb = 1e-5 * arma::speye(nb_, nb_);
#ifdef DEBUG
        data.Qb.print();
#endif
}
void ModelGenerator::construct_Qxy_lower() {
    // Qub0  --> create here, all zeros
    arma::sp_mat Qub0(nb_, ns_ * nt_);
    Qub0.zeros();
    // assemble Q.x = block_diagonal(Q.u, Q.b)
    n = size(data.Qu)[1] + size(data.Qb)[1];
    // std::cout << "n : " << n << std::endl;
    arma::sp_mat Qx(n, n);
    Qx(0, 0, size(data.Qu)) = data.Qu;
    Qx(0, nu_, size(Qub0.t())) = Qub0.t();
    Qx(nu_, 0, size(Qub0)) = Qub0;
    Qx(nu_, nu_, size(data.Qb)) = data.Qb;
    // arma::mat(Qx).submat(0,0,10,10).print();
    // Q.x|y = Q.x + t(A.x), Q.e*A.x
    arma::sp_mat Qxy(size(Qx));
    Qxy = Qx + exp(theta_[0]) * data.Ax.t() * data.Ax;
    // arma::mat(Qx).submat(0,0,10,10).print();
    data.Qxy_lower = arma::trimatl(Qxy);
    nnz = data.Qxy_lower.n_nonzero;
}

void ModelGenerator::construct_b() {
    data.b = exp(theta_[0]) * data.Ax.t() * data.y;
#ifdef DEBUG
    data.b.print();
#endif
}

void ModelGenerator::assemble_triplet_format(){
    triplets.row_idx = new size_t[nnz];
    triplets.col_ptr = new size_t[n+1];
    triplets.val = new double[nnz];
    for (int i = 0; i < nnz; ++i)
        triplets.row_idx[i] = data.Qxy_lower.row_indices[i];
    for (int i = 0; i < n + 1; ++i)
        triplets.col_ptr[i] = data.Qxy_lower.col_ptrs[i];
    for (int i = 0; i < nnz; ++i)
        triplets.val[i] = data.Qxy_lower.values[i];
}
///////////////////////////////////////////////////////////////////////////////
//                              Getter & Setter                              //
///////////////////////////////////////////////////////////////////////////////
size_t ModelGenerator::get_n(){
    return n;
}
size_t ModelGenerator::get_nnz(){
    return nnz;
}
///////////////////////////////////////////////////////////////////////////////
//                               Data Fetching                               //
///////////////////////////////////////////////////////////////////////////////
void ModelGenerator::getData() {
    // spatial matrices
    getTemporalData();
    if (nt_ > 1)
        getSpatialData();
};
void ModelGenerator::getTemporalData() {
    std::string c0_file = base_path_ + "/c0_" + ns_s + ".dat";
    std::string g1_file = base_path_ + "/g1_" + ns_s + ".dat";
    std::string g2_file = base_path_ + "/g2_" + ns_s + ".dat";
    std::string Ax_file = base_path_ + "/Ax_" + no_s + "_" + std::to_string(nu_ + nb_) + ".dat";
    std::string y_file = base_path_ + "/y_" + no_s + "_1.dat";
    if_not_exists_abort({c0_file, g1_file, g2_file, Ax_file, y_file});
    data.c0 = readCSC_sym(c0_file);
    data.g1 = readCSC_sym(g1_file);
    data.g2 = readCSC_sym(g2_file);
    data.Ax = readCSC(Ax_file);
    data.y = read_matrix(y_file, no_, 1);
};
void ModelGenerator::getSpatialData() {
    std::string g3_file = base_path_ + "/g3_" + ns_s + ".dat";
    std::string M0_file = base_path_ + "/M0_" + nt_s + ".dat";
    std::string M1_file = base_path_ + "/M1_" + nt_s + ".dat";
    std::string M2_file = base_path_ + "/M2_" + nt_s + ".dat";
    if_not_exists_abort({g3_file, M0_file, M1_file, M2_file});
    data.g3 = readCSC_sym(g3_file);
    data.M0 = readCSC_sym(M0_file);
    data.M1 = readCSC_sym(M1_file);
    data.M2 = readCSC_sym(M2_file);
#ifdef DEBUG
    arma::mat(data.g3).submat(0, 0, 10, 10).print();
    arma::mat(data.M0).submat(0, 0, nt - 1, nt - 1).print();
    arma::mat(data.M1).submat(0, 0, nt - 1, nt - 1).print();
    arma::mat(data.M2).submat(0, 0, nt - 1, nt - 1).print();
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
              << " couldn\'t be opened (not existing or failed to "
                 "open)\n";
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
    int n, nnz;
    fin >> n;
    fin >> n;
    fin >> nnz;
    arma::uvec row_ind(nnz);
    arma::uvec col_ptr(n + 1);
    arma::vec a(nnz);

    for (int i = 0; i < nnz; i++)
        fin >> row_ind[i];
    for (int i = 0; i < n + 1; i++)
        fin >> col_ptr[i];
    for (int i = 0; i < nnz; i++)
        fin >> a[i];
    fin.close();
    arma::sp_mat A_lower(row_ind, col_ptr, a, n, n);
    // create entire matrix
    arma::sp_mat A = arma::symmatl(A_lower);

    return A;
}