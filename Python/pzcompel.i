// SWIG interface
%module TPZCompEl
%include "pzsave.i"
%{
#include "../Mesh/pzcompel.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
%ignore operator<<(std::ostream &s,TPZCompEl &el);
%ignore operator << (std::ostream &out,const TPZCompElSide &celside);

%include "../Mesh/pzcompel.h"