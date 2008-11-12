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
  TPZCompMesh * cmesh = Unidimensional(1, deltaT);
  TPZGeoMesh *gmesh = cmesh->Reference();
  TPZBlackOilAnalysis an(cmesh,deltaT);

  TPZFStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ELU);
  an.SetSolver(step);

  an.SetConvergence(10, 1e-6, true);
  an.SetNewtonConvergence(100,1e-6);
  an.SetTime(deltaT);

  an.Run(cout,true);

  ofstream sol("solutionvec.txt");
  an.Solution().Print("solvec", sol);


  delete cmesh;
  delete gmesh;
  return 0;

}
