// SWIG interface
%module TPZReadGIDGrid
%{

#include "../Pre/TPZReadGIDGrid.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created

// Include the header file
//%include "pzgmesh.i"
%include "std_string.i"
%include "std_vector.i"
%include "pzstack.i"

%include "../Pre/TPZReadGIDGrid.h"

// Implement some methods need by Python
%extend MaterialDataV {
   double GETfProperties(size_t i) {
	return self->fProperties[i];
   }
   REAL __getitem__(size_t i) {
	//std::cout << "__getitem MaterialDataV__: " << std::endl;
	return self->fProperties[i];
   }
   void __setitem__(size_t i, REAL value) {
	//std::cout << "__setitem MaterialDataV__: " << std::endl;
	self->fProperties[i] = value;
   }
}








