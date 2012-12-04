// SWIG interface
%module TPZBndCond
%include "pzdiscgal.i"
%include "pzsave.i"
%{
#include "../Material/pzbndcond.h"
%}

// Ignore some methods not supported on Python and others
//%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

%include "../Material/pzbndcond.h"