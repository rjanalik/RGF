/**
 * @file      utilities.h
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

#endif
