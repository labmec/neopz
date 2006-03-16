//$Id: tetgenmesh.cpp,v 1.1 2006-03-16 13:25:52 tiago Exp $
#include "tetgenmesh.h"
#include "pzgmesh.h"
#include "pzreadtetgen.h"
#include "pzcmesh.h"
#include "pzelast3d.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzbndcond.h"
#include "pzcompel.h"
#include "TPZTimer.h"

#include "pzpoisson3d.h"

TPZCompMesh * CreateComputeMeshFromTetGen(int h, int p){

  TPZCompEl::gOrder = p;

  TPZReadTetGen read;
  TPZGeoMesh * gmesh = read.Process("malha10.1.node", "malha10.1.face", "malha10.1.ele");
  
  TPZTimer dividetime;
  dividetime.start();
  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
    int n = gmesh->NElements();
    TPZTimer dividepart;
    dividepart.start();
    for(int j = 0; j < n; j++){
      if (gmesh->ElementVec()[j]->Dimension() == 3) gmesh->ElementVec()[j]->Divide(filhos);
    }
    dividepart.stop();
    std::cout << "Divide "<< i << " - " << dividepart << std::endl;
  }
  dividetime.stop();
  std::cout << "All Divide: " << dividetime << std::endl;

  const REAL STRESS = 100.;
  const REAL EYoung = 205000.;
  const REAL Poisson= 0.3;

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);
  
#define ELASTICITY
#ifdef ELASTICITY
  TPZManVector<REAL,3> NullForce(3); 
  NullForce.Fill(0.);
  TPZElasticity3D * mat = new TPZElasticity3D(9, EYoung, Poisson, NullForce);

  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  val2.Zero();
  TPZBndCond * bcD = mat->CreateBC(0, 0,val1,val2);
  
  val2(0,0) = STRESS;
  TPZBndCond * bcN = mat->CreateBC(1, 1,val1,val2);
#endif

//#define POISSON
#ifdef POISSON
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(9, 3);
  mat->SetSymmetric();

  TPZManVector<REAL,3> convdir(3,0.);
  REAL beta = 0.0;
  mat->SetParameters(1.0, beta, convdir);
  int nstate = 1;  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bcD, *bcN;
  val2.Zero();
  bcD = mat->CreateBC(0, 0,val1,val2);
  val2(0,0) = 100.;
  bcN = mat->CreateBC(1, 1, val1, val2);
#endif
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcD);
  cmesh->InsertMaterialObject(bcN);
  
  cmesh->SetAllCreateFunctionsContinuous();
  
//   std::cout << "AutoBuild" << std::endl;
//   std::cout.flush();
  cmesh->AutoBuild();
//   std::cout << "AdjustBoundaryElements" << std::endl;
//   std::cout.flush();  
  cmesh->AdjustBoundaryElements();
//   std::cout << "CleanUpUnconnectedNodes" << std::endl;
//   std::cout.flush();  
//  cmesh->CleanUpUnconnectedNodes();
//   std::cout << "return cmesh" << std::endl;
//   std::cout.flush();  

  return cmesh;

}//CreateComputeMesh
