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
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzblock.h"
#include "TPZFrontStructMatrix.h"
#include "TPZFrontNonSym.h"

void error(char *) {}

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
// makes nothing
   result.Resize(4 /*(2 dimensões + rho + energy)*/);
   result.Fill(0.);
}

int main()
{
   //Creating the computational and geometric meshes.
   TPZFlowCompMesh * cmesh = RSCompMesh();
//   TPZFlowCompMesh * cmesh = OneElCompMesh();


// Creating the analysis object

   ofstream anFile("analysis.out");
   TPZEulerAnalysis An(cmesh, anFile);

   // Creating the structural matrix
   TPZFStructMatrix StrMatrix(cmesh);
//   TPZFrontStructMatrix<TPZFrontNonSym> StrMatrix(cmesh);
   //TPZBandStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   An.SetNewtonCriteria(1e-9, 10);
   An.SetTimeIntCriteria(1e-8,100);

   // Creating the solver for the linearized systems
   TPZStepSolver Solver;
   Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
   An.SetSolver(Solver);
//   Solver.SetGMRES(100,5,);
/*
   TPZStepSolver Pre;
   Pre.SetSSOR (2,//numiterations,
                 1.1,//overrelax
		 .0001,//tol
		 0);//From Current

   Solver.SetCG (10,//numiterations
                Pre,//pre
		1e-8,//tol
		0);//From Current*/
   An.SetSolver(Solver);

   An.Run(cout /*anFile*/);

   return 0;
}
