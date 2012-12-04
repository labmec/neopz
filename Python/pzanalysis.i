// SWIG interface
%module TPZAnalysis
%{
#include "../Analysis/pzanalysis.h"
%}

// Ignore some methods not supported on Python and others
%ignore *::Resequence(int firstel = -1);

%include "../Analysis/pzanalysis.h"