
#include "oneElementMesh.cpp"
#include "reflectedShock.cpp"
#include "reflectedShock_nonalignedmesh.cpp"
#include "SimpleShock.cpp"
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
#include "TPBSpStructMatrix.h"


void error(char *) {}

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
// makes nothing
   result.Resize(4 /*(2 dimensões + rho + energy)*/);
   result.Fill(0.);
}

int main()
{
   //Creating the computational and geometric meshes.

//   TPZFlowCompMesh * cmesh = SSCompMesh();
//   TPZFlowCompMesh * cmesh = RSCompMesh();
//   TPZFlowCompMesh * cmesh = OneElCompMesh();
   TPZFlowCompMesh * cmesh = RSNACompMesh();

// Creating the analysis object

   ofstream anFile("analysis.out");
   TPZEulerAnalysis An(cmesh, anFile);


//Solver attributes


   { // LU
   TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-10, 0);
   An.SetNewtonCriteria(1e-10, 1);
   An.SetTimeIntCriteria(1e-10,100);

   TPZStepSolver Solver;
   Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   }


/*
   { // GMRES
   TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-10, 100);
   An.SetNewtonCriteria(1e-9, 1);
   An.SetTimeIntCriteria(1e-8,10);

   //Preconditioner
   TPZStepSolver Pre;
   Pre.SetJacobi(1,//numiterations,
		 1e-9,//tol
		  0);//From Current
   Pre.SetMatrix(mat);

   //Main Solver
   TPZStepSolver Solver;
   Solver.SetGMRES(100,
		10,
		Pre,
		1e-9,
		0);
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   }
//*/
/*
   { // SSOR
   TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-10, 1000);
   An.SetNewtonCriteria(1e-9, 10);
   An.SetTimeIntCriteria(1e-8,100);

   //Main Solver
   TPZStepSolver Solver;
   Solver.SetSSOR(1, 1.1,
		1e-10,
		0);
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   }
//*/
   An.Run(cout);

   return 0;
}