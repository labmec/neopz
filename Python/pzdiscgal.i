// SWIG interface
%module TPZDiscontinuousGalerkin
%include "pzsave.i"
%{
#include "../Material/pzdiscgal.h"
%}

// Ignore some methods not supported on Python and others
//%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

%include "../Material/pzdiscgal.h"