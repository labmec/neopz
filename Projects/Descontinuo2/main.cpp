#include "oneElementMesh.cc"
#include "reflectedShock.cc"
#include "pzeuleranalysis.h"
#include "pzconslaw.h"
#include "pzmaterial.h"
#include "pzeulerconslaw.h"
#include "pzartdiff.h"
#include "pzreal.h"
#include "pzvec.h"
#include "pzflowcmesh.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include <iostream>
#include <fstream>
#include "TPZGeoElement.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzblock.h"

void error(char *) {}

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
// makes nothing
   result.Resize(4 /*(2 dimensões + rho + energy)*/);
   result.Fill(0.);
}

int main()
{
   //Creating the computational and geometric meshes.
   //TPZFlowCompMesh * cmesh = RSCompMesh();
   TPZFlowCompMesh * cmesh = OneElCompMesh();

// Creating the analysis object

   ofstream anFile("analysis.out");
   TPZEulerAnalysis An(cmesh, anFile);

   // Creating the structural matrix
   TPZBandStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   // Creating the solver for the linearized systems
   TPZStepSolver Solver;
   Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
   An.SetSolver(Solver);

   An.Run(cout /*anFile*/);

   return 0;
}
