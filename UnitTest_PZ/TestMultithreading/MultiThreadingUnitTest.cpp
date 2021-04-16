/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include <iostream>
#include <chrono>

#include "pzlog.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "pzbstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZExactFunction.h"
#include "TPZMatLaplacian.h"
#include "TPZGeoMeshTools.h"
#include "pzbndcond.h"

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz #define BOOST_TEST_MAIN pz multithreading_tests tests


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testmultithread");
#endif

#include "boost/test/unit_test.hpp"

struct SuiteInitializer
  {
    SuiteInitializer()
    {
      boost::unit_test::unit_test_log.set_threshold_level( boost::unit_test::log_warnings );
    }
};


namespace threadTest {
  constexpr int dim{2};//aux variable
  //aux function for creating 2d gmesh on unit square
  TPZGeoMesh *CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC);
  //aux function for creating 2d cmesh with laplacian mat
  TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC);

  //test the stiffness matrices in serial and parallel computations
  template <class TSTMAT>
  void CompareStiffnessMatrices(const int nThreads);

  template <class TSTMAT>
  void ComparePostProcError(const int nThreads);

  static void ForcingFunction(const TPZVec <REAL> &pt, TPZVec <STATE> &result);
  static void ExactSolution(const TPZVec <REAL> &pt, TPZVec <STATE> &sol, TPZFMatrix <STATE> &solDx);
}

BOOST_FIXTURE_TEST_SUITE(multithread_tests,SuiteInitializer)


BOOST_AUTO_TEST_CASE(multithread_assemble_test)
{
  threadTest::CompareStiffnessMatrices<TPZSkylineStructMatrix>(4);
  threadTest::CompareStiffnessMatrices<TPZBlockDiagonalStructMatrix>(4);
  threadTest::CompareStiffnessMatrices<TPZBandStructMatrix>(4);
  threadTest::CompareStiffnessMatrices<TPZSpStructMatrix>(4);
}

BOOST_AUTO_TEST_CASE(multithread_postprocerror_test)
{
  threadTest::ComparePostProcError<TPZSkylineStructMatrix>(4);
}

BOOST_AUTO_TEST_SUITE_END()


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
  auto *laplacianMat = new TPZMatLaplacian(matIdVol, threadTest::dim);

  laplacianMat->SetForcingFunction(ForcingFunction, 10);

  TPZAutoPointer<TPZFunction<STATE> > solexata = new TPZDummyFunction<STATE>(ExactSolution,10);
  TPZAutoPointer<TPZFunction<STATE>> sol(solexata);
  laplacianMat->SetExactSol(solexata);

  TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
  int bctype = 0;
  val2.Zero();
  TPZBndCond *bc = laplacianMat->CreateBC(laplacianMat, matIdBC, bctype, val1, val2);

  cmesh->InsertMaterialObject(laplacianMat);
  cmesh->InsertMaterialObject(bc);

  cmesh->SetDefaultOrder(pOrder);
  cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  return cmesh;
}

template <class TSTMAT>
void threadTest::CompareStiffnessMatrices(const int nThreads)
{
  constexpr int nDiv{4};
  constexpr int pOrder{3};
  int matIdVol;
  int matIdBC;
  auto *gMesh = CreateGMesh(nDiv, matIdVol, matIdBC);
  auto *cMesh = CreateCMesh(gMesh, pOrder, matIdVol, matIdBC);
    
  constexpr bool optimizeBandwidth{false};
  //lambda for obtaining the FE matrix
  auto GetMatrix = [cMesh, optimizeBandwidth](const int nThreads){
    TPZAnalysis an(cMesh, optimizeBandwidth);
    TSTMAT matskl(cMesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    an.Assemble();
    return an.MatrixSolver<STATE>().Matrix();
  };

  auto start = std::chrono::system_clock::now();
  auto matSerial = GetMatrix(0);
  std::cout<<typeid(TSTMAT).name()<<std::endl;
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedSerial = end - start;
  std::cout << "\tSerial time: " << elapsedSerial.count() << "s\n";

  start = std::chrono::system_clock::now();
  auto matParallel = GetMatrix(nThreads);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedParallel = end - start;
  std::cout << "\tParallel time: " << elapsedParallel.count() << "s\n";
  std::cout << "\tSpeedup: " << elapsedSerial / elapsedParallel << "x\n";

  const int nr = matParallel->Rows();
  const int nc = matParallel->Cols();
  
  TPZFMatrix<STATE> matDiff(nr,nc, 0.0);
  matSerial->Substract(matParallel, matDiff);
  const auto normDiff = Norm(matDiff);
  std::cout.precision(17);
  std::cout << std::fixed << "\tNorm diff: " << normDiff << std::endl;
  const bool checkMatNorm = IsZero(normDiff);
  BOOST_CHECK_MESSAGE(checkMatNorm,"failed");
  delete gMesh;
}

template <class TSTMAT>
void threadTest::ComparePostProcError(const int nThreads) {
  constexpr int nDiv{4};
  constexpr int pOrder{3};

  int matIdVol;
  int matIdBC;
  auto *gMesh = CreateGMesh(nDiv, matIdVol, matIdBC);
  auto *cMesh = CreateCMesh(gMesh, pOrder, matIdVol, matIdBC);

  TPZAutoPointer<TPZFunction<STATE>> solPtr(
      new TPZExactFunction<STATE>(threadTest::ExactSolution,pOrder));
  for(auto imat : cMesh->MaterialVec()){
    imat.second->SetExactSol(solPtr);
  }
  constexpr bool optimizeBandwidth{false};

  auto GetErrorVec = [cMesh, optimizeBandwidth](const int nThreads) {
    TPZAnalysis an(cMesh, optimizeBandwidth);
    TSTMAT matskl(cMesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);

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
  BOOST_CHECK_MESSAGE(pass, "failed");

  delete gMesh;
}

static void threadTest::ForcingFunction(const TPZVec <REAL> &pt, TPZVec <STATE> &result) {
  const auto x = pt[0], y = pt[1];
  result[0] = 8 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
}

static void threadTest::ExactSolution(const TPZVec <REAL> &pt, TPZVec <STATE> &sol, TPZFMatrix <STATE> &solDx) {
  const auto x = pt[0], y = pt[1];
  sol[0] = std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
  solDx(0, 0) = 2 * M_PI * std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y);
  solDx(1, 0) = 2 * M_PI * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);
  solDx(2, 0) = 0;
}
