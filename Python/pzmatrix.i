// SWIG interface
%module TPZMatrix
%{
#include "../Matrix/pzmatrix.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
%ignore *::PrintMath(const char *name = 0, std::ostream &out = std::cout);
%ignore *::AddKel;
%ignore *::AddFel;
%ignore ddot;
%ignore operator >> (std::istream& in,TPZMatrix<TT>& A) ;
//%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
//%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
//%ignore operator T*() const;

// Include the header file
%include "pzsave.i"
%include "pzvec.i"
%include "../Matrix/pzmatrix.h"

// Implement some methods need by Python
%extend TPZMatrix {
   char *__str__() {
       static char tmp[1024];
       //sprintf(tmp,"Vector(%g,%g,%g)", $self->x,$self->y,$self->z);
       sprintf(tmp,"Implement me to see the content of this matrix");
       return tmp;
   }
}

// Initializate template
%template(IntTPZMatrix) TPZMatrix<int>;
%template(DoubleTPZMatrix) TPZMatrix<double>;


