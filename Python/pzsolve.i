// SWIG interface
%module TPZSolver
%include "pzsave.i"
%{
#include "../Matrix/pzsolve.h"
%}

%module TPZMatrixSolver
%include "pzsave.i"
%{
#include "../Matrix/pzsolve.h"
%}

// Ignore some methods not supported on Python and others
//%ignore *::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

%include "../Matrix/pzsolve.h"

// Initializate template
%template(TPZSolver_REAL) TPZSolver<REAL>;
%template(TPZMatrixSolver_REAL) TPZMatrixSolver<REAL>;
//%template(TPZSolver_STATE) TPZSolver<STATE>;
