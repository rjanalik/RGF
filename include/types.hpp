/**
 * @file      types.hpp
 * @brief     Includes common type defintions included by CWC_utility.hpp, RGF.hpp, Blas.hpp
 * @date      Mon Jun 21 10:15:33 2021
 * @author    Radim
 *
 */

#ifndef __TYPES
#define __TYPES

#include <complex>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef ROWBLOCK
#define ROWBLOCK 100
#endif

#ifndef COLBLOCK
#define COLBLOCK 100
#endif

#ifndef MAX_COMM
#define MAX_COMM 20
#endif

#ifndef INF
#define INF 1E8
#endif

#ifndef SEP_TASK
#define SEP_TASK 50
#endif

#ifndef tollim
#define tollim 1E-10
#endif

#ifndef scale_norm
#define scale_norm 1.0
#endif

#ifndef __MINMAX
#define __MINMAX
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifdef Add_
#define fortran_name(x, y) (x##_)
#endif

#ifdef NoChange
#define fortran_name(x, y) (x)
#endif

#ifdef UpCase
#define fortran_name(x, y) (y)
#endif

typedef std::complex<double> CPLXType;
typedef CPLXType CPX;
typedef CPX *CPXp;
typedef CPX (*CPXpfn)(CPX);

#endif
