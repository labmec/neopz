/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include <iostream>
#include <chrono>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBuildSBFem.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSpStructMatrix.h"
#include "TPZExactFunction.h"
#include "Elasticity/TPZElasticity3D.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZGeoMeshTools.h"
#include "TPZBndCond.h"
#include "TPZSSpStructMatrix.h"
#include "TPZAnalyticSolution.h"

#include <catch2/catch.hpp>


namespace SBFemTest {

  void SBFemElasticity3D(const int nThreads);
  void SBFemBubblesDarcy(const int nThreads);
  
  TPZAutoPointer<TPZGeoMesh> CreateGMesh(const int nDiv, const int dim);
  TPZCompMesh * CreateCMesh(TPZAutoPointer<TPZGeoMesh> & gmesh, bool darcyproblem);
  void InsertMaterialElasticity3D(TPZCompMesh* cmesh);
  void InsertMaterialDarcy(TPZCompMesh* cmesh);
  void Analysis(TPZLinearAnalysis & an, const int nThreads, TPZManVector<REAL> & errorVec);
  
  // Material Data
  constexpr int EGroup{15};
  constexpr int ESkeleton{12};
  constexpr int EMat1{1};
  constexpr int EBc1{-1};

  // Solution data
  constexpr int pOrder{3};
  TLaplaceExample1 LaplaceExact;
  TElasticity3DAnalytic ElastExact;
}

TEST_CASE("SBFEM convergence test","[sbfem][analysis]")
{
  std::cout << "#######################\nTesting SBFEM-Elasticity 3D approximation\n";
// SBFemTest::SBFemElasticity3D(4);
  std::cout << "\n\n#######################\nTesting SBFEM-Bubbles approximation\n";
  SBFemTest::SBFemBubblesDarcy(32);
}

TPZAutoPointer<TPZGeoMesh> SBFemTest::CreateGMesh(const int nDiv, const int dim)
{
  MMeshType meshType;
  if(dim == 2){
    meshType = MMeshType::EQuadrilateral;
  } else{
    meshType = MMeshType::EHexahedral;
  }
  
  static TPZManVector<REAL,3> minX(3,-1);
  static TPZManVector<REAL,3> maxX(3,1);
  
  TPZVec<int> nDivs(dim,nDiv);

  TPZManVector<int,5> matIdVec(1+dim*2, SBFemTest::EBc1);
  matIdVec[0] = SBFemTest::EGroup;

  return TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,true);
}

TPZCompMesh * SBFemTest::CreateCMesh(TPZAutoPointer<TPZGeoMesh> & gmesh, bool darcyproblem)
{
  auto *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);

  std::map<int, int> matmap;
  matmap[SBFemTest::EGroup] = SBFemTest::EMat1;

  TPZBuildSBFem build(gmesh, SBFemTest::ESkeleton, matmap);
  build.StandardConfiguration();
  build.DivideSkeleton(0);

  if(darcyproblem){
    InsertMaterialDarcy(cmesh);
  } else{
    InsertMaterialElasticity3D(cmesh);
  }
  
  build.BuildComputationMesh(*cmesh);

  return cmesh;
}

void SBFemTest::InsertMaterialDarcy(TPZCompMesh * cmesh)
{
  TPZFMatrix<STATE> val1(1, 1, 0.);
  const TPZManVector<double> val2(1, 0.);
  auto *mat = new TPZDarcyFlow(SBFemTest::EMat1, 2);
  auto forcingFunction = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
    LaplaceExact.ForcingFunction()->Execute(x, u);
  };
  mat->SetForcingFunction(forcingFunction, pOrder);
  mat->SetConstantPermeability(1.);
  cmesh->InsertMaterialObject(mat);

  auto bc = mat->CreateBC(mat, SBFemTest::EBc1, 0, val1, val2);
  cmesh->InsertMaterialObject(bc);
  bc->SetForcingFunctionBC(LaplaceExact.ExactSolution());

  auto bcs = mat->CreateBC(mat, SBFemTest::ESkeleton, 1, val1, val2);
  cmesh->InsertMaterialObject(bcs);
}

void SBFemTest::InsertMaterialElasticity3D(TPZCompMesh * cmesh)
{
  TPZFMatrix<STATE> val1(3, 3, 0.);
  const TPZManVector<double> val2(3, 0.);
  auto *mat = new TPZElasticity3D(SBFemTest::EMat1);
  mat->SetBigNumber(1.e20);
  mat->SetMaterialDataHook(ElastExact.fE, ElastExact.fPoisson);
  cmesh->InsertMaterialObject(mat);

  auto bc = mat->CreateBC(mat, SBFemTest::EBc1, 0, val1, val2);
  bc->SetForcingFunctionBC(ElastExact.ExactSolution());
  cmesh->InsertMaterialObject(bc);

  auto bcs = mat->CreateBC(mat, SBFemTest::ESkeleton, 1, val1, val2);
  cmesh->InsertMaterialObject(bcs);
}

void SBFemTest::Analysis(TPZLinearAnalysis & an, const int nThreads, TPZManVector<REAL> &errorVec)
{
  TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matssp(an.Mesh());
  matssp.SetNumThreads(nThreads);
  an.SetStructuralMatrix(matssp);

  TPZStepSolver <STATE> *direct = new TPZStepSolver<STATE>;
  direct->SetDirect(ELDLt);
  an.SetSolver(*direct);
  delete direct;
  direct = 0;

  an.Assemble();
  an.Solve();
  an.SetThreadsForError(nThreads);

  an.PostProcessError(errorVec, false);
};

void SBFemTest::SBFemElasticity3D(const int nThreads){
  TPZSBFemElementGroup::gDefaultPolynomialOrder = 0;
  ElastExact.fProblemType = TElasticity3DAnalytic::ELoadedBeam;
  ElastExact.fE = 1.; ElastExact.fPoisson = 0.3;

  constexpr int nDiv{1};
  constexpr int dim{3};

  TPZAutoPointer<TPZGeoMesh> gMesh = CreateGMesh(nDiv, dim);
  auto *cmesh = CreateCMesh(gMesh, false);

  TPZLinearAnalysis an(cmesh, true);
  an.SetExact(ElastExact.ExactSolution());

  TPZManVector<REAL> errorVecPar;
  auto start = std::chrono::system_clock::now();
  Analysis(an, nThreads, errorVecPar);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedParallel = end - start;
  std::cout << "\tParallel time: " << elapsedParallel.count() << "s\n";
  std::cout << "\tParallel errors: " << errorVecPar << std::endl;

  bool pass = true;
  for (auto i = 0; i < 3; i++) {
    if (fabs(errorVecPar[i]) > 1e-10) {
      pass = false;
    }
  }
  REQUIRE(pass);

  delete cmesh;
}

void SBFemTest::SBFemBubblesDarcy(const int nThreads) {
  TPZSBFemElementGroup::gDefaultPolynomialOrder = 3;
  constexpr int nDiv{4};
  constexpr int pOrder{3};
  constexpr int dim{2};
  LaplaceExact.fExact = TLaplaceExample1::ECosCos;

  TPZAutoPointer<TPZGeoMesh> gMesh = CreateGMesh(nDiv, dim);
  auto *cmesh = CreateCMesh(gMesh, true);

  TPZLinearAnalysis an(cmesh, true);
  an.SetExact(SBFemTest::LaplaceExact.ExactSolution());

  auto start = std::chrono::system_clock::now();
  TPZManVector<REAL> errorVecSer;
  Analysis(an, 0, errorVecSer);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedSerial = end - start;
  std::cout << "\tSerial time: " << elapsedSerial.count() << "s\n";
  std::cout << "\tSerial errors: " << errorVecSer << std::endl;

  start = std::chrono::system_clock::now();
  TPZManVector<REAL> errorVecPar;
  Analysis(an, nThreads, errorVecPar);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedParallel = end - start;
  std::cout << "\tParallel time: " << elapsedParallel.count() << "s\n";
  std::cout << "\tParallel errors: " << errorVecPar << std::endl;

  std::cout << "\tSpeedup: " << elapsedSerial / elapsedParallel << "x\n";

  bool pass = true;
  for (auto i = 0; i < errorVecSer.size(); i++) {
    if (!IsZero(errorVecSer[i] - errorVecPar[i])) {
      pass = false;
    }
  }
  REQUIRE(pass);

  delete cmesh;
}
