// SWIG interface
%module TPZManVector
%{
#include "../Util/pzmanvector.h"
#include "../Pre/TPZReadGIDGrid.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
//%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
//%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
//%ignore operator T*() const;

// Include the header file
%include "../Util/pzmanvector.h"
%include "pzvec.i"

// Implement some methods need by Python
%extend TPZManVector {
   char *__str__() {
       static char tmp[1024];
       //sprintf(tmp,"Vector(%g,%g,%g)", $self->x,$self->y,$self->z);
       sprintf(tmp,"Implement me to see the content of this TPZManVector");
       return tmp;
   }
   T __getitem__(size_t i) {
std::cout << "__getitem TPZManVector__: " << std::endl;
	return (*self)[ i ];
   }
   void __setitem__(size_t i, T value) {
	(*self)[ i ] = value;
   }
   TPZManVector <T> copy() {
	TPZManVector <T> tmp (*self);
	return tmp;
   }
}

// Initializate template
//%template(IntTPZManVector) TPZManVector<int>;
//%template(DoubleTPZManVector) TPZManVector<double>;
%template(MaterialDataVTPZManVector) TPZManVector< MaterialDataV,200 >;
%template(REALTPZManVector) TPZManVector< REAL,200 >;
%include "std_string.i"
%template(TPZManVector_string) TPZManVector< std::string, 10 >;


