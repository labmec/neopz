#ifndef TPZLapackH
#define TPZLapackH

#ifdef USING_LAPACK
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
typedef __CLPK_doublecomplex vardoublecomplex;
typedef __CLPK_complex varfloatcomplex;

#elif USING_MKL
#include <mkl.h>
typedef MKL_Complex16 vardoublecomplex;
typedef MKL_Complex8 varfloatcomplex;


#elif WIN32
#include <complex>

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

typedef lapack_complex_double vardoublecomplex;
typedef lapack_complex_float varfloatcomplex;

#include "cblas.h"
#include "lapacke.h"

#else
#include <clapack.h>

#endif
#endif

#endif
