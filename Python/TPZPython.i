// SWIG interface
%module TPZPython
%{
typedef double REAL; //This is the default configuration
typedef double STATE; //This is the default configuration
%}

typedef double REAL; //This is the default configuration
typedef double STATE; //This is the default configuration

%include "pzvec.i"
%include "pzfilebuffer.i"
%include "pzsave.i"
%include "pzstack.i"
%include "pzfmatrix.i"
%include "pzmatrix.i"
%include "pzstrmatrix.i"
%include "pzcompel.i"
%include "pzgeoel.i"
%include "pzgmesh.i"
%include "pzcmesh.i"
%include "pzpoisson3d.i"
%include "pzmaterial.i"
%include "pzdiscgal.i"
%include "tpzreadgidgrid.i"
%include "pzmanvector.i"
%include "pzanalysis.i"
%include "pzskylstrmatrix.i"
%include "pzstepsolver.i"
%include "pzsolve.i"
%include "pzlog.i"
%include "pzbndcond.i"
//%include "log4cxx.i"
%include "pzchunk.i"
%include "pzadmchunk.i"
