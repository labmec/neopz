// SWIG interface
%module TPZGeoEl
%include "pzsave.i"
%{
#include "../Mesh/pzgeoel.h"
%}

// Ignore some methods not supported on Python and others
//%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
%ignore TPZGeoEl::ComputePermutationNormals(int, TPZVec<int>&);

%include "../Mesh/pzgeoel.h"