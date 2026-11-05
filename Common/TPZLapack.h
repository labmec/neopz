#ifndef TPZLapackH
#define TPZLapackH

#ifdef USING_LAPACK
#ifdef MKLLAPACK
#include <mkl.h>
typedef MKL_Complex16 vardoublecomplex;
typedef MKL_Complex8 varfloatcomplex;
//lapack_int already defined typedef long long lapack_int;

#elif MACOSX
#include <Accelerate/Accelerate.h>
typedef __CLPK_doublecomplex vardoublecomplex;
typedef __CLPK_complex varfloatcomplex;
typedef int lapack_int;
#else
#include <complex>

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

typedef lapack_complex_double vardoublecomplex;
typedef lapack_complex_float varfloatcomplex;
typedef int lapack_int;
#include "lapacke.h"
#endif
#endif

#endif
