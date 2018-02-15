// SWIG interface
%module TPZSavable
%{
#include "../Save/pzsave.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::Restore(TPZStream&, void*);
%ignore *::Register(int, TPZSavable* (*)(TPZStream&, void*));
%ignore *::Compare(TPZSavable *,bool);
%ignore *::Compare(TPZSavable *);
%ignore *::Compare(TPZSavable *,int);
%ignore *::Write(TPZStream &,int);
%ignore *::Write(int const *);
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
//%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
//%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
//%ignore operator T*() const;

// Include the header file
%include "pzfilebuffer.i"
%include "../Save/pzsave.h"
%include "pzvec.i"

// Implement some methods need by Python
%extend TPZSavable {
   char *__str__() {
       static char tmp[1024];
       //sprintf(tmp,"Vector(%g,%g,%g)", $self->x,$self->y,$self->z);
       sprintf(tmp,"Implement me to see the content of this matrix");
       return tmp;
   }
}

// Initializate template
//%template(IntTPZSrtMatrix) TPZSavable<int>;
//%template(DoubleTPZSrtMatrix) TPZSavable<double>;


