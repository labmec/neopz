#ifndef TPZLapackH
#define TPZLapackH

#ifdef USING_LAPACK
#ifdef MKLLAPACK
#include <mkl.h>
typedef MKL_Complex16 vardoublecomplex;
typedef MKL_Complex8 varfloatcomplex;

#elif MACOSX
#include <Accelerate/Accelerate.h>
typedef __CLPK_doublecomplex vardoublecomplex;
typedef __CLPK_complex varfloatcomplex;

#else
#include <complex>

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

typedef lapack_complex_double vardoublecomplex;
typedef lapack_complex_float varfloatcomplex;

#include "lapacke.h"
#endif
#endif

#endif
