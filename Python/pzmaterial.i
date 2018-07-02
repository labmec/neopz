// SWIG interface
%module TPZMaterial
%include "pzsave.i"
%{
#include "../Material/TPZMaterial.h"
%}

// Ignore some methods not supported on Python and others
//%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

%include "../Material/pzmaterial.h"