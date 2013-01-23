// SWIG interface
%module TPZFMatrix
%{
#include "../Matrix/pzfmatrix.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
//%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
//%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
//%ignore operator T*() const;
%ignore *::operator-() const;
%ignore operator const double *;

// Include the header file
%include "pzsave.i"
%include "pzmatrix.i"
%include "../Matrix/pzfmatrix.h"

// Implement some methods need by Python
%extend TPZFMatrix {
   char *__str__() {
       static char tmp[1024];
       //sprintf(tmp,"Vector(%g,%g,%g)", $self->x,$self->y,$self->z);
       sprintf(tmp,"Implement me to see the content of this matrix");
       return tmp;
   }
}

// Initializate template
//%template(IntTPZFMatrix) TPZFMatrix<int>;
//%template(DoubleTPZFMatrix) TPZFMatrix<double>;
%template(TPZFMatrix_REAL) TPZFMatrix<REAL>;


