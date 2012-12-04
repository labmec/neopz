// SWIG interface
%module TPZMatPoisson3d
%include "pzsave.i"
%include "pzmaterial.i"
%include "pzdiscgal.i"
%{
#include "../Material/pzpoisson3d.h"
%}

// Ignore some methods not supported on Python and others
//%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

%include "../Material/pzpoisson3d.h"

// Implement some methods need by Python
%extend TPZMatPoisson3d {
}