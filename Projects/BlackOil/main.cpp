//Classes utilit�ias
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "malha.h"
#include "pzblackoilanalysis.h"

using namespace std;

int main(){

  const double deltaT = 1.;
  TPZCompMesh * cmesh = Unidimensional(0, deltaT);
  TPZGeoMesh *gmesh = cmesh->Reference();
  TPZBlackOilAnalysis an(cmesh,deltaT);

  TPZFStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ELU);
  an.SetSolver(step);

  an.SetConvergence(1000, 1e-6, true);
  an.SetNewtonConvergence(100,1e-6);
  an.SetTime(deltaT);

  TPZFMatrix initialSol = an.Solution();
  for(int i = 0; i < initialSol.Rows()/2; i++){
    initialSol(2*i,0)   = 26.e6;
    initialSol(2*i+1,0) = 1.;
  }
  an.SetInitialSolution(initialSol);

  an.Run(cout,false);

  ofstream sol("solutionvec.txt");
  an.Solution().Print("solvec", sol);


  delete cmesh;
  delete gmesh;
  return 0;

}
