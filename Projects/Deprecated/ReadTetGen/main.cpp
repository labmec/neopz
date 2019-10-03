//$Id: main.cpp,v 1.3 2007-06-08 00:05:31 cesar Exp $

/**
 * Validation test of TPZReadtetGen class wich implements the interface between tetgen and PZ.
 * March 03, 2006
 */

#include "pzgmesh.h"
#include "pzreadtetgen.h"
#include "pzcmesh.h"
#include "pzelast3d.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzbndcond.h"
#include "pzcompel.h"


TPZCompMesh * CreateComputeMesh(TPZGeoMesh & gmesh);

int main(){

  TPZReadTetGen read;
  TPZGeoMesh * gmesh = read.Process("malha10.1.node", "malha10.1.face", "malha10.1.ele");
//  gmesh->Print(std::cout);

  TPZCompEl::SetgOrder(1);
  TPZCompMesh * cmesh = CreateComputeMesh(*gmesh);
  TPZAnalysis an(cmesh);

  /*TPZFrontStructMatrix <TPZFrontSym> */TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ECholesky);
  an.SetSolver(step);
  std::cout << "\nNEquations = " << cmesh->NEquations();  std::cout.flush();
  std::cout << "\nAssemble..."; std::cout.flush();
  an.Assemble();
  std::cout << "\nSolve..."; std::cout.flush();
  an.Solve();
  std::cout << "\nFinished"; std::cout.flush();

}//main

TPZCompMesh * CreateComputeMesh(TPZGeoMesh & gmesh){

  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.0;

  TPZCompMesh *cmesh = new TPZCompMesh(&gmesh);
  cmesh->SetDimModel(3);

  TPZManVector<REAL,3> NullForce(3);
  NullForce.Fill(0.);
  TPZElasticity3D * mat = new TPZElasticity3D(9, EYoung, Poisson, NullForce);

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  TPZBndCond * bcD = mat->CreateBC(0, 0,val1,val2);

  val2(0,0) = STRESS;
  TPZBndCond * bcN = mat->CreateBC(1, 1,val1,val2);

  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);

  cmesh->SetAllCreateFunctionsContinuous();

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  return cmesh;

}//CreateComputeMesh
