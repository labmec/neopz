// SWIG interface
%module TPZSkylineStructMatrix
%{
#include "../StrMatrix/pzskylstrmatrix.h"
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
%include "pzstrmatrix.i"
%include "../StrMatrix/pzskylstrmatrix.h"

// Implement some methods need by Python


// Initializate template
//%template(IntTPZFMatrix) TPZFMatrix<int>;
//%template(DoubleTPZFMatrix) TPZFMatrix<double>;
//%template(TPZFMatrix_REAL) TPZFMatrix<REAL>;


