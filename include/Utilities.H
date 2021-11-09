/**
 * @file      utilities.hpp
 * @brief     TODO
 * \todo add brief description of utlities.h
 * @date      Mon Jun 21 10:22:38 2021
 * @author    Radim
 * @copyright
 *
 * This module
 */

#ifndef __UTILITIES
#define __UTILITIES

#include "Blas.H"
#include <cmath>
#include <ctime>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif


void icopy(int, int *, int *);

int Round(double);

double get_time(double);

/************************************************************************************************/

template <typename T> void init_var(T *var, int N);

template <typename T> void init_var(T *var, int N) {
  for (int i = 0; i < N; i++) {
    var[i] = (T)0;
  }
}

/************************************************************************************************/

template <typename T, typename W> T convert(W val);

template <> inline double convert(CPX val) { return real(val); }

template <> inline CPX convert(CPX val) { return val; }

template <> inline double convert(double val) { return val; }

template <> inline CPX convert(double val) { return CPX(val, 0.0); }

/************************************************************************************************/
void parse_args(int argc, char *argv[], std::string &base_path, size_t &ns,
                size_t &nt, size_t &nb, size_t &no,
                std::ostream &stream = std::cout);

namespace utilities {

  void print_header(std::string title = "", size_t length = 100, char symbol = '=', size_t sep_width = 5, char sep = ' ', std::ostream &stream = std::cout);
  size_t mf_block_index(size_t supernode, size_t *diag_pos, size_t b_size);
  size_t mf_block_lda(size_t supernode, size_t ns, size_t nt, size_t nd);
  void print_matrix_structure(T *M, size_t m, size_t n, bool RMO = true);
  void print_matrix_structure_with_fixed_effects(T *M, size_t m, size_t n, size_t nd, bool RMO = true);
  void print_MF_structure(T *MF, size_t ns, size_t nt, size_t nd, size_t *diag_pos);
  size_t number_of_digits(double n);
  void get_max_digits_per_column(const T *M, size_t n, size_t *max_len_per_col);
  void print_matrix(T *M,  size_t m, size_t n, bool RMO = true);
  void print_csr(size_t *ia, size_t *ja, double *a, size_t n, size_t nd, bool structureOnly = false);
  void print_ia_ja_a(const size_t *ia, const size_t *ja, const T *a, size_t n);
}

#endif
