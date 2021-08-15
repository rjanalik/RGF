#ifndef __MODELGENERATOR_H_
#define __MODELGENERATOR_H_

#include <armadillo>
#include <string>
#include <vector>
#include "utilities.hpp"

class ModelGenerator {
public:
  ModelGenerator();
  ModelGenerator(size_t ns, size_t nt, size_t nb, size_t no, arma::vec theta,
                 std::string base_path);
  void getData();
  bool file_exists(const std::string &file_name);
  size_t get_n();
  size_t get_nnz();
  void construct_model();
  void assemble_triplet_format();
  struct Triplet {
    size_t *ja;
    size_t *ia;
    double *a;
    ~Triplet() {
      delete[] ja;
      delete[] ia;
      delete[] a;
    }
  } triplets;

private:
  void getTemporalData();
  void getSpatialData();
  void construct_Qu();
  void construct_Qu_spatial(arma::vec theta_u);
  void construct_Qu_spatio_temporal(arma::vec theta_u);
  void construct_Qb();
  void construct_Qx();
  void construct_Qxy_lower();
  void construct_b();
  void if_not_exists_abort(std::string const file_name);
  void if_not_exists_abort(std::initializer_list<std::string> const file_names);
  arma::sp_mat readCSC(std::string filename);
  arma::sp_mat readCSC_sym(std::string const filename);
  arma::mat read_matrix(std::string filename, int n_row, int n_col);

  size_t ns_;
  size_t nt_;
  size_t nb_;
  size_t no_;
  size_t nu_;
  size_t n_;
  size_t nnz_;
  arma::vec theta_;
  std::string base_path_;
  struct Data {
    // spatial matrices;
    arma::sp_mat c0;
    arma::sp_mat g1;
    arma::sp_mat g2;
    arma::sp_mat Ax;
    arma::vec y;
    // temporal matrices;
    arma::sp_mat g3;
    arma::sp_mat M0;
    arma::sp_mat M1;
    arma::sp_mat M2;
    // model matrices
    arma::sp_mat Qu;
    arma::sp_mat Qb;
    arma::sp_mat Qxy_lower;
    arma::vec b;

  } data;
};
#endif // __MODELGENERATOR_H_
