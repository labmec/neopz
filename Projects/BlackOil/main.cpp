//Classes utilit�ias
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sstream>

#include "malha.h"
#include "pzblackoilanalysis.h"

using namespace std;

int main(){

/** Dados da simulacao */
  /** Passo de tempo */
  const double deltaT = 2629743.8*1.;///1 month
  /** Aresta do dominio em planta */
  const double Lr = 1524.;
  /** Altura do reservatorio */
  const double Hr = 9.144;
  /** Raio do poco */
  const double rw = 0.3;
  /** Numero de divisoes otimo ao redor dos pocos */
  const double pi = 4.*atan(1.);
  const int ndivOpt = log(4.*(Lr-rw)/(pi*rw))/log(2.)+0.5;
  /** Numero de divisoes adotado ao redor dos pocos */
  const int ndiv = 14-11;//+6+3;
  cout << "\nNDiv opt = " << ndivOpt << "  ndiv adotado = " << ndiv << "\n";
  /** Numero de divisoes na espessura */
  const int nlayers = 2;

//   TPZCompMesh * cmesh = Unidimensional(0, deltaT);   TPZGeoMesh *gmesh = cmesh->Reference();
//   TPZCompMesh * cmesh = UnidimensionalGravidade(1, deltaT);
 TPZGeoMesh * gmesh = QuarterFiveSpot(ndiv, nlayers, Lr,  Hr, 2.*rw, 10.);
//   TPZGeoMesh * gmesh = QuarterFiveSpotReg(Lr, Hr);
  DivideMalha(gmesh); /*DivideMalha(gmesh);*/
//   DivideTornoPocos(gmesh);DivideTornoPocos(gmesh);
  TPZCompMesh * cmesh = QuarterFiveSpot(gmesh, deltaT);

  cout << "\nNEquations = " << cmesh->NEquations() << "\n"; cout.flush();
  TPZBlackOilAnalysis an(cmesh,deltaT);

  /*TPZFStructMatrix*/TPZBandStructMatrix/*TPZFrontStructMatrix<TPZFrontNonSym>*/ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ELU);
  an.SetSolver(step);

  an.SetConvergence(20.*12.*2629743.8/deltaT+0.5, 1e-6, true);
  an.SetNewtonConvergence(1000,1e-5);
  an.SetSaveFrequency(1,0);

  TPZVec<string> scalnames(3);
  TPZVec<string> vecnames(0);
//   scalnames[0] = "WaterPressure";
//   scalnames[1] = "OilPressure";
//   scalnames[2] = "WaterSaturation";
//   scalnames[3] = "OilSaturation";
  scalnames[0] = "OilPressure";
  scalnames[1] = "WaterSaturation";
  scalnames[2] = "OilSaturation";

  std::stringstream filedx;
  filedx << "result.dx";
  an.DefineGraphMesh(3, scalnames, vecnames, &(filedx.str()[0]));
  an.SetTime(deltaT);

  TPZFMatrix initialSol = an.Solution();
  for(int i = 0; i < initialSol.Rows()/2; i++){
    initialSol(2*i,0)   = 33094840.;
    initialSol(2*i+1,0) = 1.;
  }

//    initialSol(0,0) = 30.e6;
//    initialSol(1,0) = 1.;
//    initialSol(2,0) = 21.e6;
//    initialSol(3,0) = 1.;

  an.SetInitialSolution(initialSol);

  an.Run(cout,/*false*/true);

  ofstream sol("solutionvec.txt");
  an.Solution().Print("solvec", sol);


  delete cmesh;
  delete gmesh;
  return 0;

}
