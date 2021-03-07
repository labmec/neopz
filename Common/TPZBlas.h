#ifndef TPZBlasH
#define TPZBlasH

#ifdef USING_BLAS
#ifdef MKLBLAS
#include <mkl.h>
#else
#include "cblas.h"
#endif
#endif

#endif
