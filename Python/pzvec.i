// SWIG interface
%module TPZVec
%{
#include "../Util/pzvec.h"
#include "../Pre/TPZReadGIDGrid.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
%ignore operator T*() const;

// Include the header file
%include "../Util/pzvec.h"

// Implement some methods need by Python
%extend TPZVec {
   char *__str__() {
       static char tmp[1024];
       //sprintf(tmp,"Vector(%g,%g,%g)", $self->x,$self->y,$self->z);
       sprintf(tmp,"Implement me to see the content of this TPZVec");
       return tmp;
   }
   T __getitem__(size_t i) {
//std::cout << "__getitem TPZVec__: " << std::endl;
	return (*self)[ i ];
   }
   void __setitem__(size_t i, T value) {
	(*self)[ i ] = value;
   }
   TPZVec <T> copy() {
	TPZVec <T> tmp (*self);
	return tmp;
   }
}

// Initializate template
//%template(IntTPZVec) TPZVec<int>;
//%template(DoubleTPZVec) TPZVec<double>;
%include "pzpoisson3d.i"
%include "pzmaterial.i"
%template(TPZVec_TPZMatPoisson3d) TPZVec<TPZMatPoisson3d *>;
%template(TPZVec_TPZMaterial) TPZVec<TPZMaterial *>;
%template(TPZVec_REAL) TPZVec<REAL>;
%template(TPZVec_MaterialDataV) TPZVec< MaterialDataV >;
%include "std_string.i"
%template(TPZVec_string) TPZVec< std::string >;
%include "pzgeoel.i"
%template(TPZVec_TPZGeoEl) TPZVec<TPZGeoEl *>;

