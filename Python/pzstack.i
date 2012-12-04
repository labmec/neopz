// SWIG interface
%module TPZStack
%{
#include "../Util/pzstack.h"
#include "../Pre/TPZReadGIDGrid.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::operator=; //method copy created
%ignore *::operator[]; //methods __getitem__ and __setitem__ created
%ignore operator<<(std::ostream &out, const std::pair<int,int> &element);
%ignore operator<<( std::ostream& Out, const TPZVec< T2 >& v );
%ignore operator T*() const;

// Include the header file
%include "pzsave.i"
%include "pzvec.i"
%include "pzmanvector.i"
%include "../Util/pzstack.h"
%include "../Pre/TPZReadGIDGrid.h"

// Implement some methods need by Python
%extend TPZStack {
   T __getitem__(size_t i) {
//std::cout<< "AQUI - getitem [" << i << "]" << std::endl;
	return (*self)[ i ];
   }
   void __setitem__(size_t i, T value) {
	(*self)[ i ] = value;
//std::cout<< "AQUI - setitem [" << i << "] " << std::endl;
   }
   T GETA(size_t i) {
	return (*self)[ i ];
   }
}


// Initializate template
%include "std_vector.i"
%template(REALVector) std::vector<REAL>;
//%template(IntTPZStack) TPZStack<int>;
//%template(DoubleTPZStack) TPZStack<double>;
%template(REALTPZStack) TPZStack<REAL>;
%template(MatDataTPZStack) TPZStack<MaterialDataV>;


