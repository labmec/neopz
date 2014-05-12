//
//  CohesiveTests.cpp
//  PZ
//
//  Created by Nathan Shauer on 9/5/14.
//  Copyright (c) 2014 LabMec-Unicamp. All rights reserved.
//

#include <iostream>

#include "CohesiveTests.h"


#include "pzcmesh.h"
#include "pzanalysis.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"

//Teste CohesiveBC
#include "pznlelasmat.h"
#include "pznonlinanalysis.h"
#include "TPZCohesiveBC.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
//Teste CohesiveBC


void GetSolAtLeft(TPZCompMesh	*cmesh)
{
	TPZGeoMesh *gmesh = cmesh->Reference();
	TPZGeoEl *gel = gmesh->ElementVec()[0];
	if (!gel ) {
		DebugStop();
	}
	TPZGeoElSide gelside(gel,0);
	TPZGeoElSide gelsideneig = gelside.Neighbour();
	TPZCompEl *cel = NULL;
	while (gelsideneig != gelside) {
		TPZGeoEl *gelRefinado = gelsideneig.Element();
		cel = gelRefinado->Reference();
		if (cel) {
			break;
		}
		gelsideneig = gelsideneig.Neighbour();
	}

	TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *> (cel);
	if (!sp) {
		DebugStop();
	}
	int var = 9;
	TPZManVector<REAL,3> Solout(3,0.), qsi(3,-1.);
	sp->Solution(qsi, var, Solout);
	std::cout << "uy = " << Solout[1] << std::endl;
}

void ElastNLTestWithCohesive()
{
  TPZGeoMesh *gmesh = CreateGeoMeshCohe();
  TPZCompMesh *cmesh = CreateCMeshCohe(gmesh);
  TPZNonLinearAnalysis an(cmesh,std::cout);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  TPZSkylineStructMatrix skyl(cmesh);
  an.SetStructuralMatrix(skyl);
  an.SetSolver(step);
  
  SolveNLElasticity(cmesh,an);
	GetSolAtLeft(cmesh);
	

  int dim = 2;
  TPZStack<std::string> scalnames,vecnames;
  vecnames.push_back("Strain");
  vecnames.push_back("Displacement");
  scalnames.push_back("SigmaX");
  scalnames.push_back("SigmaY");
  
  an.DefineGraphMesh(dim, scalnames, vecnames, "ElastNLSol.vtk");
  
  an.PostProcess(0);
  
  int nsteps = 1;
  for (int i = 2; i <= nsteps ; i++) {
   /*
    it = cmesh->MaterialVec().begin();
    sigma = i;
    for (; it != cmesh->MaterialVec().end(); it++) {
      TPZCohesiveBC *bc = dynamic_cast<TPZCohesiveBC*> (it->second);
      if (bc) {
        bc->SetCohesiveData(sigma, sigma, sigma);
      }
    }
    */
    SolveNLElasticity(cmesh,an);
    an.PostProcess(0);
  }
  
}

void SolveNLElasticity(TPZCompMesh *cmesh, TPZNonLinearAnalysis &an)
{
 
  REAL tol = 1.e-9;
  int numiter = 10;
  an.IterativeProcess(std::cout, tol, numiter);

}


TPZGeoMesh* CreateGeoMeshCohe()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  int nnodes = 4, nodeid = 0;
  
  TPZVec<REAL> coord(3,0.);
  gmesh->NodeVec().Resize(nnodes);
  
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
  
  coord[0] = 1.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
  
  coord[1] = 1.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
  
  coord[0] = 0.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
  
  TPZVec<long> TopolQuad(4,0);
  TopolQuad[0] = 0;
  TopolQuad[1] = 1;
  TopolQuad[2] = 2;
  TopolQuad[3] = 3;
  
  long index = 0;
  TPZGeoEl *gel = NULL;
	int matid = 1;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matid, index);
  gel->SetId(index);
  index++;
  
//  long bcdir = -1;
//  gel->CreateBCGeoEl(1,bcdir);
  
  long bcfixedx = -3;
  gel->CreateBCGeoEl(2, bcfixedx);
  
  long bcneu = -2;
  gel->CreateBCGeoEl(6,bcneu);

//  long bcneu2 = -4;
//  gel->CreateBCGeoEl(4,bcneu2);
  
  TPZVec<long> TopolLine(2,0);
  TopolLine[0] = 0;
  TopolLine[1] = 1;
  int cohebc = 2;
  gel = gmesh->CreateGeoElement(EOned, TopolLine, cohebc, index);
  gel->SetId(index);
  index++;
  
  gmesh->BuildConnectivity();
  
  int nref = 3;
  TPZVec<TPZGeoEl *> sons;
  for (int iref = 0; iref < nref; iref++) {
    int nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
      TPZGeoEl *gel = gmesh->ElementVec()[iel];
      if (!gel->HasSubElement()) {
        gel->Divide(sons);
      }
    }
  }
  
  std::ofstream out("ElastGmeshCohe.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
  
  return gmesh;
}

TPZCompMesh* CreateCMeshCohe(TPZGeoMesh *gmesh)
{
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  /// criar materiais
	int dim = 2;
	
  TPZVec<REAL> force(dim,0.);
  
  //int planestress = 1;
  int planestress = 1;
  
  REAL ela = 40., nu = 0., fx = 0., fy = 0.;
	int matid = 1;
  TPZNLElasticityMaterial * material1 = new TPZNLElasticityMaterial(matid,
                                                                    ela,
                                                                    nu,
                                                                    fx,
                                                                    fy,
                                                                    planestress);
  
  TPZMaterial * mat1(material1);
  
  ///criar malha computacional
  cmesh->SetDefaultOrder(2);
	cmesh->SetDimModel(dim);
  cmesh->InsertMaterialObject(mat1);
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  int bcdir = -1, bcneu = -2, bcfixedx = -3;
  int dir = 0, neu = 1, fixedx = 3;
//  TPZMaterial * BCond11 = material1->CreateBC(mat1, bcdir, dir, val1, val2);
//  cmesh->InsertMaterialObject(BCond11);

  val2(1,0) = 1;
  TPZMaterial * BCond12 = material1->CreateBC(mat1, bcneu, neu, val1, val2);
  cmesh->InsertMaterialObject(BCond12);
  
  val2.Zero();
//  TPZMaterial * BCond13 = material1->CreateBC(mat1, bcfixedx, fixedx, val1, val2);
//  cmesh->InsertMaterialObject(BCond13);
	
//	val1.Zero();
//  val2.Zero();
//  val2(1,0) = -0.00001;
//	int bcneu2 = -4;
//  TPZMaterial * BCond14 = material1->CreateBC(mat1, bcneu2, neu, val1, val2);
//  cmesh->InsertMaterialObject(BCond14);
  
  int cohesiveid = 2;
  TPZCohesiveBC * material2 = new TPZCohesiveBC(cohesiveid);
  const REAL SigmaT = 2., DeltaC = 1, DeltaT = 0.2;
  material2->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  TPZMaterial *mat2(material2);
	cmesh->InsertMaterialObject(mat2);
  
  /*
  val1.Redim(2,2);
  val2.Redim(2,1);
  val2(1,0) = 1.;
  TPZMaterial * BCond21 = material1->CreateBC(mat1, bcneu, neu, val1, val2);
	cmesh->InsertMaterialObject(BCond21);
  */
   
  //cmesh->SetAllCreateFunctionsContinuous();
  cmesh->SetAllCreateFunctionsContinuousWithMem();
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

void ElastTest()
{
  TPZGeoMesh *gmesh = CreateGeoMesh();
  TPZCompMesh *cmesh = CreateCMesh(gmesh);
	for (int i = 0; i < 300; i++) {
	  SolveLinearElasticity(cmesh);
	}
  //SolveLinearElasticity(cmesh);
	
	delete cmesh;
	delete gmesh;
}

void ElastNLTest()
{
  TPZGeoMesh *gmesh = CreateGeoMesh();
  TPZCompMesh *cmesh = CreateCMesh(gmesh);
  TPZNonLinearAnalysis an(cmesh,std::cout);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  TPZSkylineStructMatrix skyl(cmesh);
  an.SetStructuralMatrix(skyl);
  an.SetSolver(step);
  
  REAL sigma = 1;
  
  SolveNLElasticity(cmesh,an);
  int nsteps = 1;
	
  int dim = 2;
  TPZStack<std::string> scalnames,vecnames;
  vecnames.push_back("Strain");
  vecnames.push_back("Displacement");
  scalnames.push_back("SigmaX");
  scalnames.push_back("SigmaY");
  
  an.DefineGraphMesh(dim, scalnames, vecnames, "ElastNLSol.vtk");
	
	an.PostProcess(0);
	
	delete cmesh;
	delete gmesh;
}


void SolveLinearElasticity(TPZCompMesh *cmesh)
{
  
  TPZAnalysis an(cmesh);
  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);
  TPZParFrontStructMatrix<TPZFrontSym<STATE> > skyl(cmesh);
	//TPZSkylineStructMatrix skyl(cmesh);
	skyl.SetNumThreads(6);
	
	an.SetStructuralMatrix(skyl);
  an.SetSolver(step);
  an.Run();
  /*
  int dim = 2;
  TPZStack<std::string> scalnames,vecnames;
  vecnames.push_back("Strain");
  vecnames.push_back("Displacement");
  scalnames.push_back("SigmaX");
  scalnames.push_back("SigmaY");
  
  an.DefineGraphMesh(dim, scalnames, vecnames, "ElastSol.vtk");
  an.PostProcess(0);
	 */
	
}


TPZGeoMesh* CreateGeoMesh()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  int nnodes = 4, nodeid = 0;
  
  TPZVec<REAL> coord(3,0.);
  gmesh->NodeVec().Resize(nnodes);
  
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
  
  coord[0] = 1.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;

  coord[1] = 2.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;

  coord[0] = 0.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
  
  TPZVec<long> TopolQuad(4,0);
  TopolQuad[0] = 0;
  TopolQuad[1] = 1;
  TopolQuad[2] = 2;
  TopolQuad[3] = 3;
  
  long index = 0;
  TPZGeoEl *gel = NULL;
	int matid = 1;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matid, index);
  gel->SetId(index);
  index++;

  long bcdir = -1, bcneu = -2;
  
  gel->CreateBCGeoEl(4,bcdir);
  gel->CreateBCGeoEl(6,bcneu);
  
  gmesh->BuildConnectivity();
  
  int nref = 2;
  TPZVec<TPZGeoEl *> sons;
  for (int iref = 0; iref < nref; iref++) {
    int nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
      TPZGeoEl *gel = gmesh->ElementVec()[iel];
      if (!gel->HasSubElement()) {
        gel->Divide(sons);
      }
    }
  }
  
  std::ofstream out("ElastGmesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
  
  return gmesh;
}

TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh)
{
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  /// criar materiais
	int dim = 2;
	
  TPZVec<REAL> force(dim,0.);
  
  //int planestress = 1;
  int planestrain = 1;
  
  REAL ela = 40., nu = 0., fx = 0., fy = 0.;
	int matid = 1;
  TPZNLElasticityMaterial * material1 = new TPZNLElasticityMaterial(matid,
                                                                ela,
                                                                nu,
                                                                fx,
                                                                fy,
                                                                planestrain);

  TPZMaterial * mat1(material1);
  
  ///criar malha computacional
  cmesh->SetDefaultOrder(2);
	cmesh->SetDimModel(dim);
  cmesh->InsertMaterialObject(mat1);
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  int bcdir = -1, bcneu = -2;
  int dir = 0, neu = 1;  
	TPZMaterial * BCond11 = material1->CreateBC(mat1, bcdir, dir, val1, val2);
  cmesh->InsertMaterialObject(BCond11);
  
  val1.Redim(2,2);
  val2.Redim(2,1);
  val2(1,0) = 1.;
  TPZMaterial * BCond21 = material1->CreateBC(mat1, bcneu, neu, val1, val2);
  
  cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(BCond21);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}