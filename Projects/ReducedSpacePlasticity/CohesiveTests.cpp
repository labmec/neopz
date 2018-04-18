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
  var = sp->Material()->VariableIndex("SigmaY");
  sp->Solution(qsi, var, Solout);
  std::cout << "SigmaY = " << Solout[0] << std::endl;
}


void ElastNLTestWithCohesive()
{
  TPZGeoMesh *gmesh = CreateGeoMeshCohe();
  REAL val = 0.3;
  TPZCompMesh *cmesh = CreateCMeshCohe(gmesh,val);
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
}

void CohesiveTwoLoads()
{
  TPZGeoMesh *gmesh = CreateGeoMeshCohe();
  REAL val = 0.3;
  TPZCompMesh *cmesh = CreateCMeshCohe(gmesh,val);
  TPZNonLinearAnalysis an(cmesh,std::cout);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  TPZSkylineStructMatrix skyl(cmesh);
  an.SetStructuralMatrix(skyl);
  an.SetSolver(step);
  
  SolveNLElasticity(cmesh,an);
	GetSolAtLeft(cmesh);
	
  // setando para atualizar o material coesivo na primeira iteracao do newton
  std::map<int,TPZMaterial*> MatMap = cmesh->MaterialVec();
  std::map<int,TPZMaterial*>::iterator it = MatMap.find(2); // id do material coesivo
  if (it == MatMap.end()) {
    DebugStop(); // Material coesivo nao esta na malha?
  }
  else{
    TPZCohesiveBC *matcohe = dynamic_cast<TPZCohesiveBC*>(it->second);
    matcohe->SetUpdateMem(true);
  }
  it = MatMap.find(-4); // id of the recalque
  if (it == MatMap.end()) {
    DebugStop(); // Material coesivo nao esta na malha?
  }
  else{
    TPZBndCond *bnd = dynamic_cast<TPZBndCond*>(it->second);
    TPZFMatrix<> val2(2,1,0.);
    val2(1,0) = 0.25;
    bnd->Val2() = val2;
  }
  
  
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
}

void SolveNLElasticity(TPZCompMesh *cmesh, TPZNonLinearAnalysis &an)
{
	int iter = 0;
	REAL error = 1.e10;
	int numeq = cmesh->NEquations();
	
	TPZFMatrix<STATE> prevsol(an.Solution());
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
	
  REAL tol = 5.e-4;
  int numiter = 20;
	while(error > tol && iter < numiter) {
		
		an.Solution().Redim(0,0);
		an.Assemble();
    
    //assim que fiz o primeiro assemble atualizando, nao devo mais atualizar o maerial coesivo
    // POSSO FAZER O PRIMEIRO ASSEMBLE FORA ao inves de seta toda santa vez (fica pra quando for fazer o verdadeiro)
    std::map<int,TPZMaterial*> MatMap = cmesh->MaterialVec();
    std::map<int,TPZMaterial*>::iterator it = MatMap.find(2); // id do material coesivo
    if (it == MatMap.end()) {
      DebugStop(); // Material coesivo nao esta na malha?
    }
    else{
      TPZCohesiveBC *matcohe = dynamic_cast<TPZCohesiveBC*>(it->second);
      matcohe->SetUpdateMem(false);
    }
    
		an.Solve();
    an.Solution() += prevsol;
    
    
		prevsol -= an.Solution();
		REAL normDeltaSol = Norm(prevsol);
		prevsol = an.Solution();
		an.LoadSolution(prevsol);
		an.AssembleResidual();
		double NormResLambda = Norm(an.Rhs());
		double norm = NormResLambda;
    std::cout << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << std::endl;
		
		if(norm < tol) {
      std::cout << "\nTolerancia atingida na iteracao : " << (iter+1) << std::endl;
      std::cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << std::endl << std::endl;
			
		} else
			if( (norm - error) > 1.e-9 ) {
        std::cout << "\nDivergent Method\n";
			}
		error = norm;
		iter++;
    std::cout.flush();
	}
 
  /*
  REAL tol = 1.e-8;
  int numiter = 30;
  an.IterativeProcess(std::cout, tol, numiter);
   */

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
  
  TPZVec<int64_t> TopolQuad(4,0);
  TopolQuad[0] = 0;
  TopolQuad[1] = 1;
  TopolQuad[2] = 2;
  TopolQuad[3] = 3;
  
  int64_t index = 0;
  TPZGeoEl *gel = NULL;
	int matid = 1;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matid, index);
  gel->SetId(index);
  index++;
  
//  int64_t bcdir = -1;
//  gel->CreateBCGeoEl(1,bcdir);
  
  int64_t bcfixedx = -3;
  gel->CreateBCGeoEl(2, bcfixedx);
  
  int64_t bcneu = -2;
  gel->CreateBCGeoEl(6,bcneu);
  
  int64_t bcrecalque = -4;
  gel->CreateBCGeoEl(6,bcrecalque);
  
  TPZVec<int64_t> TopolLine(2,0);
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

TPZCompMesh* CreateCMeshCohe(TPZGeoMesh *gmesh, REAL val)
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
  int dir = 0, neu = 1, fixedx = 3, bcrecalque = -4;;
  
//  TPZMaterial * BCond11 = material1->CreateBC(mat1, bcdir, dir, val1, val2);
//  cmesh->InsertMaterialObject(BCond11);

//  val2(1,0) = 1;
//  TPZMaterial * BCond12 = material1->CreateBC(mat1, bcneu, neu, val1, val2);
//  cmesh->InsertMaterialObject(BCond12);
  
  val2.Zero();
//  TPZMaterial * BCond13 = material1->CreateBC(mat1, bcfixedx, fixedx, val1, val2);
//  cmesh->InsertMaterialObject(BCond13);

  val2(1,0) = val;
  TPZMaterial * BCond14 = material1->CreateBC(mat1, bcrecalque, dir, val1, val2);
  cmesh->InsertMaterialObject(BCond14);
  
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
  
  TPZVec<int64_t> TopolQuad(4,0);
  TopolQuad[0] = 0;
  TopolQuad[1] = 1;
  TopolQuad[2] = 2;
  TopolQuad[3] = 3;
  
  int64_t index = 0;
  TPZGeoEl *gel = NULL;
	int matid = 1;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matid, index);
  gel->SetId(index);
  index++;

  int64_t bcdir = -1, bcneu = -2;
  
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

void OpenFractureTest()
{
  TPZGeoMesh *gmesh = CreateGeoMeshToOpen();
  REAL val = 0.3;
  TPZCompMesh *cmesh = CreateCMeshToOpen(gmesh,val);
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
  
  an.DefineGraphMesh(dim, scalnames, vecnames, "Propagatio.vtk");
  
  an.PostProcess(0);
	
  UpdateBcValue(cmesh,0.55);  // setando para atualizar o material coesivo na primeira iteracao do newton
  SolveNLElasticity(cmesh,an);
  an.PostProcess(0);
  
  UpdateBcValue(cmesh,0.7);  // setando para atualizar o material coesivo na primeira iteracao do newton
  SolveNLElasticity(cmesh,an);
  an.PostProcess(0);
  
 	GetSolAtLeft(cmesh);
  

}

TPZGeoMesh* CreateGeoMeshToOpen()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  int nnodes = 6, nodeid = 0;
  
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
  
  coord[0] = 2.;
  coord[1] = 0.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;

  coord[0] = 2.;
  coord[1] = 1.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;

  
  TPZVec<int64_t> TopolQuad(4,0);
  TopolQuad[0] = 0;
  TopolQuad[1] = 1;
  TopolQuad[2] = 2;
  TopolQuad[3] = 3;
  
  int64_t index = 0;
  TPZGeoEl *gel = NULL;
	int matid = 1;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matid, index);
  gel->SetId(index);
  index++;
  
  int64_t bcrecalque = -4;
  gel->CreateBCGeoEl(6,bcrecalque);

  TopolQuad[0] = 1;
  TopolQuad[1] = 4;
  TopolQuad[2] = 5;
  TopolQuad[3] = 2;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matid, index);
  gel->SetId(index);
  index++;

  
  int64_t bcdir = -1;
  gel->CreateBCGeoEl(4,bcdir);
  
  int64_t bcfixedx = -3;
  gel->CreateBCGeoEl(2, bcfixedx);
  
//  int64_t bcneu = -2;
//  gel->CreateBCGeoEl(6,bcneu);
  
  //int64_t bcrecalque = -4;
  gel->CreateBCGeoEl(6,bcrecalque);
  
  TPZVec<int64_t> TopolLine(2,0);
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
  
  std::ofstream out("ElastGmeshPropag.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
  
  return gmesh;
}

TPZCompMesh* CreateCMeshToOpen(TPZGeoMesh *gmesh, REAL val)
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
  int dir = 0, neu = 1, fixedx = 3, bcrecalque = -4;;
  
  TPZMaterial * BCond11 = material1->CreateBC(mat1, bcdir, dir, val1, val2);
  cmesh->InsertMaterialObject(BCond11);
  
  //  val2(1,0) = 1;
  //  TPZMaterial * BCond12 = material1->CreateBC(mat1, bcneu, neu, val1, val2);
  //  cmesh->InsertMaterialObject(BCond12);
  
  val2.Zero();
  //TPZMaterial * BCond13 = material1->CreateBC(mat1, bcfixedx, fixedx, val1, val2);
  //cmesh->InsertMaterialObject(BCond13);
  
  val2(1,0) = val;
  TPZMaterial * BCond14 = material1->CreateBC(mat1, bcrecalque, dir, val1, val2);
  cmesh->InsertMaterialObject(BCond14);
  
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

void UpdateBcValue(TPZCompMesh *cmesh, REAL val)
{
  std::map<int,TPZMaterial*> MatMap = cmesh->MaterialVec();
  std::map<int,TPZMaterial*>::iterator it = MatMap.find(2); // id do material coesivo
  if (it == MatMap.end()) {
    DebugStop(); // Material coesivo nao esta na malha?
  }
  else{
    TPZCohesiveBC *matcohe = dynamic_cast<TPZCohesiveBC*>(it->second);
    matcohe->SetUpdateMem(true);
  }
  it = MatMap.find(-4); // id of the recalque
  if (it == MatMap.end()) {
    DebugStop(); // Material coesivo nao esta na malha?
  }
  else{
    TPZBndCond *bnd = dynamic_cast<TPZBndCond*>(it->second);
    TPZFMatrix<> val2(2,1,0.);
    val2(1,0) = val;
    bnd->Val2() = val2;
  }
}

void TestContinuousDiscontinuous()
{
	TPZGeoMesh *gmesh = CreateGeomeshContDisc();
	TPZCompMesh *cmesh = CreateCmeshContDisc(gmesh);
	
	std::ofstream out2("gMeshContDisc.txt");
	cmesh->Reference()->Print(out2);
	
	std::ofstream out("cMeshContDisc.txt");
	cmesh->Print(out);
	
	delete cmesh;
	delete gmesh;
}

TPZGeoMesh* CreateGeomeshContDisc()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  int nnodes = 6, nodeid = 0;
  
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
  
  coord[0] = 2.;
  coord[1] = 0.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
	
  coord[0] = 2.;
  coord[1] = 1.;
  gmesh->NodeVec()[nodeid].SetNodeId(nodeid);
  gmesh->NodeVec()[nodeid].SetCoord(coord);
  nodeid++;
	
	// id das cond contorno
	int bcdircont = -1;
 	int bcdirdisc = -2;
	int matidcont = 1;
	int matiddisc = 2;
  
	// Elemento 0
  TPZVec<int64_t> TopolQuad(4,0);
  TopolQuad[0] = 0;
  TopolQuad[1] = 1;
  TopolQuad[2] = 2;
  TopolQuad[3] = 3;
  
  int64_t index = 0;
  TPZGeoEl *gel = NULL;
	int matid = 1;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matidcont, index);
  gel->SetId(index);
  index++;
  
	// Cond contorno do El 0
  gel->CreateBCGeoEl(4,bcdircont);
	gel->CreateBCGeoEl(6,bcdircont);
	gel->CreateBCGeoEl(7,bcdircont);
	
	// Elemento 1
  TopolQuad[0] = 1;
  TopolQuad[1] = 4;
  TopolQuad[2] = 5;
  TopolQuad[3] = 2;
  gel = gmesh->CreateGeoElement(EQuadrilateral, TopolQuad, matiddisc, index);
  gel->SetId(index);
  index++;

	// Cond contorno El 1
  gel->CreateBCGeoEl(4,bcdirdisc);
	gel->CreateBCGeoEl(5,bcdirdisc);
	gel->CreateBCGeoEl(6,bcdirdisc);
	
	gmesh->AddInterfaceMaterial(matidcont,matiddisc,matiddisc); //Adicionar um material de interface associados aos elementos mat2 e mat1 do material.
	gmesh->AddInterfaceMaterial(matiddisc,matidcont,matiddisc);
   
  gmesh->BuildConnectivity();
  
	//Refinamento
	/*
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
	*/
  
  std::ofstream out("ElastGmeshContDisc.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
  
  return gmesh;
}

TPZCompMesh* CreateCmeshContDisc(TPZGeoMesh *gmesh)
{
	
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
	
	/// criar materiais
	int dim = 2;
	
  TPZVec<REAL> force(dim,0.);
  
  //int planestress = 1;
  int planestress = 1;
  
  REAL ela = 40., nu = 0., fx = 0., fy = 0.;
	int bcdircont = -1;
 	int bcdirdisc = -2;
	int matidcont = 1;
	int matiddisc = 2;
  TPZNLElasticityMaterial * material1 = new TPZNLElasticityMaterial(matidcont,
                                                                    ela,
                                                                    nu,
                                                                    fx,
                                                                    fy,
                                                                    planestress);
	TPZNLElasticityMaterial * material2 = new TPZNLElasticityMaterial(matiddisc,
                                                                    ela,
                                                                    nu,
                                                                    fx,
                                                                    fy,
                                                                    planestress);
	
  
  TPZMaterial * mat1(material1);
	TPZMaterial * mat2(material2);
  
  ///criar malha computacional
  cmesh->SetDefaultOrder(2);
	cmesh->SetDimModel(dim);
	
  cmesh->InsertMaterialObject(mat1);
  cmesh->InsertMaterialObject(mat2);
	
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  int dir = 0, neu = 1;
  
  TPZMaterial * BCond11 = material1->CreateBC(mat1, bcdircont, dir, val1, val2);
  cmesh->InsertMaterialObject(BCond11);
	
	TPZMaterial * BCond12 = material2->CreateBC(mat2, bcdirdisc, dir, val1, val2);
  cmesh->InsertMaterialObject(BCond12);
  
  //  val2(1,0) = 1;
  //  TPZMaterial * BCond12 = material1->CreateBC(mat1, bcneu, neu, val1, val2);
  //  cmesh->InsertMaterialObject(BCond12);
 
	
	std::set<int> MatCont, MatDisc;
	MatCont.insert(matidcont);
	MatCont.insert(bcdircont);
	MatDisc.insert(matiddisc);
	MatDisc.insert(bcdirdisc);
	
  //cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetAllCreateFunctionsContinuous();
	
	cmesh->AutoBuild(MatCont);
	gmesh->ResetReference();

	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->AutoBuild(MatDisc);
	cmesh->LoadReferences(); 

  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	
  return cmesh;
	
}