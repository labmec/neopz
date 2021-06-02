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
#include "pzskylstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "pzbstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZExactFunction.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZGeoMeshTools.h"
#include "TPZBndCond.h"

#include <catch2/catch.hpp>


namespace threadTest {
  constexpr int dim{2};//aux variable
  //aux function for creating 2d gmesh on unit square
  TPZGeoMesh *CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC);
  //aux function for creating 2d cmesh with laplacian mat
  TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC);

  template <class TSTMAT>
  void ComparePostProcError(const int nThreads);

  constexpr int solOrder{4};
  auto ExactSolution = [](const TPZVec<REAL> &loc,
         TPZVec<STATE>&u,
         TPZFMatrix<STATE>&gradU){
        const auto &x=loc[0];
        const auto &y=loc[1];
        u[0]=(x*x-1)*(y*y-1);
        gradU(0,0) = 2*x*(y*y-1);
        gradU(1,0) = 2*y*(x*x-1);
        gradU(2,0) = 0;//optional
  };
  constexpr int rhsPOrder{2};
  const auto ForcingFunction = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
        const REAL &x = loc[0];
        const REAL &y = loc[1];
        u[0] = 2*y*y+2*x*x-4;
        u[0] *= -1;
  };
}

TEST_CASE("Parallel post process error test","[multithread_tests][multithread][struct][analysis]")
{
  threadTest::ComparePostProcError<TPZSkylineStructMatrix<STATE>>(4);
}


TPZGeoMesh *threadTest::CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC)
{
  constexpr MMeshType meshType = MMeshType::ETriangular;
  
  static TPZManVector<REAL,3> minX(3,0);
  static TPZManVector<REAL,3> maxX(3,1);
  maxX[2] = 0.;
  TPZVec<int> nDivs(dim,nDiv);
  TPZManVector<int,5> matIdVec(5, -1);
  matIdVec[0] = 1;
  matIdVol = matIdVec[0];
  matIdBC = matIdVec[1];
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(threadTest::dim,minX,maxX,matIdVec,nDivs,meshType,true);
}

TPZCompMesh *threadTest::CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC)
{
  auto *cmesh = new TPZCompMesh(gmesh);
  auto *poissonMat = new TPZMatPoisson(matIdVol, threadTest::dim);

  poissonMat->SetForcingFunction(ForcingFunction, rhsPOrder);

  TPZFNMatrix<1, REAL> val1(1, 1, 0.);
  TPZVec<STATE> val2(1, 0.);
  int bctype = 0;
  TPZBndCond *bc = poissonMat->CreateBC(poissonMat, matIdBC, bctype, val1, val2);

  cmesh->InsertMaterialObject(poissonMat);
  cmesh->InsertMaterialObject(bc);

  cmesh->SetDefaultOrder(pOrder);
  cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  return cmesh;
}

template <class TSTMAT>
void threadTest::ComparePostProcError(const int nThreads) {
  constexpr int nDiv{4};
  constexpr int pOrder{3};

  int matIdVol;
  int matIdBC;
  auto *gMesh = CreateGMesh(nDiv, matIdVol, matIdBC);
  auto *cMesh = CreateCMesh(gMesh, pOrder, matIdVol, matIdBC);
  constexpr bool optimizeBandwidth{false};
  auto GetErrorVec = [cMesh, optimizeBandwidth](const int nThreads) {
    TPZLinearAnalysis an(cMesh, optimizeBandwidth);
    TSTMAT matskl(cMesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    an.SetExact(threadTest::ExactSolution);
    TPZStepSolver <STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;

    an.Assemble();
    an.Solve();

    
    an.SetThreadsForError(nThreads);

    TPZManVector <REAL> errorVec(3, 0.);
    int64_t nelem = cMesh->NElements();
    cMesh->ExpandSolution();
    cMesh->ElementSolution().Redim(nelem, 10);

    an.PostProcessError(errorVec);
    return errorVec;
  };

  auto start = std::chrono::system_clock::now();
  auto errorVecSer = GetErrorVec(0);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedSerial = end - start;
  std::cout << "\tSerial time: " << elapsedSerial.count() << "s\n";
  std::cout << "\tSerial errors: " << errorVecSer << std::endl;

  start = std::chrono::system_clock::now();
  auto errorVecPar = GetErrorVec(nThreads);
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

  delete cMesh;
  delete gMesh;
}
