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

#include "Blas.hpp"
#include <armadillo>
#include <cmath>
#include <ctime>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

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
}


#endif