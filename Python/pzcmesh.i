// SWIG interface
%module TPZCompMesh
%include "pzsave.i"
%{
#include "../Mesh/pzcmesh.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::Consolidate();
%ignore *::ComputeMesh(TPZVec<int> &accumlist,int numaggl);
%ignore *::Check();
%ignore *::ConnectVec();
%ignore *::Block();

%include "../Mesh/pzcmesh.h"