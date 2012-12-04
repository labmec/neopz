// SWIG interface
%module TPZGeoMesh
%include "pzsave.i"
%{
#include "../Mesh/pzgmesh.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::ElementVec();
%ignore *::NodeVec();

%include "../Mesh/pzgmesh.h"