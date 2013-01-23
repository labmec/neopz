// SWIG interface
%module TPZGeoMesh
%include "pzsave.i"
%include "pzgeoel.i"
%{
#include "../Mesh/pzgmesh.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::ElementVec();
%ignore *::NodeVec();

%include "../Mesh/pzgmesh.h"