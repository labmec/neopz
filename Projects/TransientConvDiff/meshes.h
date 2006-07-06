// -*- c++ -*-

//$Id: meshes.h,v 1.1 2006-07-06 15:48:25 tiago Exp $

#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include <set>

//4 quadrados iguais, podendo ser divididos. Podem ser continuos ou descontinuos
TPZCompMesh *CreateMesh(int h);

// u = x*(x-1.)*y*(y-1.)*T*exp(-T);
void ForcingFunction(TPZVec<REAL> &pto, TPZVec<REAL> &force);
void ExactSolution(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv);

// u = x*(x-1.)*y*(y-1.)*T;
void ForcingFunction2(TPZVec<REAL> &pto, TPZVec<REAL> &force);
void ExactSolution2(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv);
