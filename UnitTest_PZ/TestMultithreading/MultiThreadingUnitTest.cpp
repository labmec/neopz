/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include <iostream>

#include "pzlog.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZMatLaplacian.h"
#include "TPZGeoMeshTools.h"

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

BOOST_FIXTURE_TEST_SUITE(multithread_tests,SuiteInitializer)
    
BOOST_AUTO_TEST_CASE(multithread_assemble_test)
{
  constexpr int dim{2};
  constexpr int nDiv{4};
  constexpr MMeshType meshType{MMeshType::ETriangular};
  constexpr int pOrder{3};
  TPZManVector<int,1> matIdVec(1);
  TPZGeoMesh *gMesh = [&]() -> TPZGeoMesh *
  {
    static TPZManVector<REAL,2> minX(2,0);
    static TPZManVector<REAL,2> maxX(2,1);
    TPZVec<int> nDivs(dim,nDiv);
    return TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,false);
  }();

  TPZCompMesh *cMesh = [&]() -> TPZCompMesh *
  {    
    auto *cmesh = new TPZCompMesh(gMesh);
 
    auto *laplacianMat = new TPZMatLaplacian(matIdVec[0], dim);
    cmesh->InsertMaterialObject(laplacianMat);
 
    cmesh->SetDefaultOrder(pOrder);
 
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
  }();
  constexpr bool optimizeBandwidth{true};
  //lambda for obtaining the matrix
  auto GetMatrix = [&](const int nThreads){
    TPZAnalysis an(cMesh, optimizeBandwidth);
    TPZSkylineStructMatrix matskl(cMesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
    return an.Solver().Matrix();
    
  };
  
  auto matSerial = GetMatrix(0);
  auto matParallel = GetMatrix(4);
  std::cout<< "r = " << matParallel->Rows() << std::endl;
  std::cout<< "c = " << matParallel->Cols() << std::endl;
  TPZFMatrix<STATE> matDiff;
  matSerial->Substract(matParallel.operator*(),matDiff);
  const auto normDiff = Norm(matDiff);
  std::cout << normDiff << std::endl;
  const bool checkMatNorm= IsZero(normDiff);
  BOOST_CHECK_MESSAGE(checkMatNorm,"failed");
  delete cMesh;
  delete gMesh;
}


BOOST_AUTO_TEST_SUITE_END()
