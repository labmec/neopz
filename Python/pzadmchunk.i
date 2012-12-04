// SWIG interface
%module TPZAdmChunkVector
%{
#include "../Util/pzadmchunk.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
//%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
//%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
//%ignore operator T*() const;

// Include the header file
%include "../Util/pzadmchunk.h"
%include "pzvec.i"
%include "pzmanvector.i"
%include "pzchunk.i"

// Implement some methods need by Python
%extend TPZAdmChunkVector {
   char *__str__() {
       static char tmp[1024];
       //sprintf(tmp,"Vector(%g,%g,%g)", $self->x,$self->y,$self->z);
       sprintf(tmp,"Implement me to see the content of this TPZAdmChunkVector");
       return tmp;
   }
   T __getitem__(size_t i) {
	return (*self)[ i ];
   }
   void __setitem__(size_t i, T value) {
	(*self)[ i ] = value;
   }
   TPZAdmChunkVector <T> copy() {
	TPZAdmChunkVector <T> tmp (*self);
	return tmp;
   }
}

%template(IntTPZAdmChunkVector) TPZAdmChunkVector<int>;
%include "pzgeoel.i"
%template(TPZAdmChunkVector_TPZGeoEl) TPZAdmChunkVector<TPZGeoEl *>;


