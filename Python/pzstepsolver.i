// SWIG interface
%module TPZStepSolver
%{
#include "../Matrix/pzstepsolver.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
//%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
//%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
//%ignore operator T*() const;

// Include the header file
%include "pzsave.i"
%include "pzmatrix.i"
%include "pzsolve.i"
%include "../Matrix/pzstepsolver.h"

// Implement some methods need by Python


// Initializate template
%template(TPZStepSolver_REAL) TPZStepSolver<REAL>;
//%template(TPZStepSolver_STATE) TPZStepSolver<STATE>;



