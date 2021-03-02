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

#include "TPZMatLaplacian.h"
#include "TPZGeoMeshTools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz #define BOOST_TEST_MAIN pz multithreading_tests tests


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testmultithread"));
#endif

#include "boost/test/unit_test.hpp"

struct SuiteInitializer
  {
    SuiteInitializer()
    {
      InitializePZLOG();
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
};

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

  TLaplaceExample1* exact = new TLaplaceExample1;
  exact->fExact = TLaplaceExample1::ESinSin;
  laplacianMat->SetForcingFunctionExact(exact->Exact());
  laplacianMat->SetForcingFunction(exact->ForcingFunction());

  TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
  int bctype = 0;
  val2.Zero();
  TPZBndCond *bc = laplacianMat->CreateBC(laplacianMat, matIdBC, bctype, val1, val2);
  bc->TPZMaterial::SetForcingFunction(exact->Exact());

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
    return an.Solver().Matrix();
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
void threadTest::ComparePostProcError(const int nThreads)
{
    constexpr int nDiv{4};
    constexpr int pOrder{3};

    int matIdVol;
    int matIdBC;
    auto *gMesh = CreateGMesh(nDiv, matIdVol, matIdBC);
    auto *cMesh = CreateCMesh(gMesh, pOrder, matIdVol, matIdBC);

    constexpr bool optimizeBandwidth{false};

    auto GetErrorVec = [cMesh, optimizeBandwidth](const int nThreads){
        TPZAnalysis an(cMesh, optimizeBandwidth);
        TSTMAT matskl(cMesh);
        matskl.SetNumThreads(nThreads);
        an.SetStructuralMatrix(matskl);

        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;

        an.Assemble();
        an.Solve();

        TLaplaceExample1* exact = new TLaplaceExample1;
        exact->fExact = TLaplaceExample1::ESinSin;
        an.SetExact(exact->ExactSolution());
        an.SetThreadsForError(nThreads);

        TPZManVector<REAL> errorVec(3, 0.);
        int64_t nelem = cMesh->NElements();
        cMesh->LoadSolution(cMesh->Solution());
        cMesh->ExpandSolution();
        cMesh->ElementSolution().Redim(nelem, 10);

        an.PostProcessError(errorVec);
        return errorVec;
    };

    auto errorVecSer = GetErrorVec(0);
    std::cout << "Serial errors: " << errorVecSer << std::endl;

    auto errorVecPar = GetErrorVec(nThreads);
    std::cout << "Parallel errors: " << errorVecPar << std::endl;

    bool pass = true;
    for (auto i = 0; i < errorVecSer.size(); i++) {
        if (!IsZero(errorVecSer[i] - errorVecPar[i])) {
            pass = false;
        }
    }
    BOOST_CHECK_MESSAGE(pass, "failed");

    delete gMesh;
}
