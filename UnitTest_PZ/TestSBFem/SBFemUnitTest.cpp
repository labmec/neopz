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
  void VerifyShape(TPZCompMesh *cmesh);
  
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
  SBFemTest::SBFemElasticity3D(4);
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

  // Verifying sbfem bubble functions
  VerifyShape(cmesh);

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


void SBFemTest::VerifyShape(TPZCompMesh *cmesh) {
  bool pass = true;
  int64_t nel = cmesh->NElements();
  TPZFMatrix<STATE> Sol(cmesh->NEquations(),1,0.);

  // change to true if you wanna display the element outputs
  bool outputs = false;

  for (size_t el = 0; el < nel; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    TPZSBFemElementGroup *celgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
    if(celgr) {
      if(outputs) {
        std::cout << "element group " << el << std::endl;
      }
      auto &phinv = celgr->PhiInverse();
      int sz = phinv.Rows();
      TPZFMatrix<STATE> phinvreal(sz,sz);
      for(int i=0; i<sz; i++) for(int j=0; j<sz; j++) phinvreal(i,j) = phinv(i,j).real();
      int nc = cel->NConnects();
      auto &elgr = celgr->GetElGroup();
      auto nelgr = elgr.size();
      for (size_t i = 0; i < nelgr; i++)
      {
        TPZCompEl *celvol = elgr[i];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(celvol);
        if(sbfem) {
          TPZManVector<REAL,3> xi(2,0.),sol(1,0.);
          //xi[1] = -1.;
          TPZFMatrix<REAL> phi,dphixi;
          sbfem->Shape(xi,phi,dphixi);
          if(outputs) {
            sbfem->Print(std::cout);
            std::cout << "values of phi \n";
            for(int i = 0; i<sz; i++) std::cout << phi(i,0) << " ";
            std::cout << std::endl;
          }
          int firstshape = 0;
          for(int ic = 0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seq = c.SequenceNumber();
            int64_t pos = cmesh->Block().Position(seq);
            int64_t blsize = cmesh->Block().Size(seq);
            for(int ishape = 0; ishape < blsize; ishape++) {
              Sol.Zero();
              Sol(pos+ishape,0) = 1.;
              cmesh->LoadSolution(Sol);
              sbfem->Solution(xi,1,sol);
              if (outputs){
                if(0 &&ic == nc-1)
                {
                  TPZVec<STATE> coefs(sz);
                  for(int i=0; i<sz; i++) coefs[i] = phinvreal(i,ishape+firstshape);
                  std::cout << " column in phinv " << coefs << std::endl;
                  auto coefload = celgr->CoeficientsReal();
                  int sz = coefload.Rows();
                  TPZVec<STATE> coefloadvec(sz);
                  for(int i=0; i<sz; i++) coefloadvec[i] = coefload(i,0);
                  std::cout << "working coeficients " << coefloadvec << std::endl;
                }
                std::cout << "subel " << i << " shape " << ishape+firstshape << " phi(ishape) " << phi(firstshape+ishape,0) << " sol " << sol[0] << std::endl;
              }
              if (!IsZero(phi(firstshape+ishape,0)-sol[0]))
              {
                pass = false;
              }
            }
            firstshape += blsize;
          }
        }
      }
    }
  }
  std::cout << "\n Testing SBFEM Bubble functions\n";
  REQUIRE(pass);
  if(pass)
    std::cout << "Bubble shape functions PASSED\n\n" << std::endl;
}