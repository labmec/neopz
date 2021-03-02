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


namespace threadTest{
  constexpr int dim{2};//aux variable
  //aux function for creating 2d gmesh on unit square
  TPZGeoMesh *CreateGMesh(const int nDiv, int& matIdVol);
  //aux function for creating 2d cmesh with laplacian mat
  TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol);


  //test the stifness matrices in serial and parallel computations
  template <class TSTMAT>
  void CompareStiffnessMatrices(const int nThreads);
};

BOOST_FIXTURE_TEST_SUITE(multithread_tests,SuiteInitializer)


BOOST_AUTO_TEST_CASE(multithread_assemble_test)
{
  threadTest::CompareStiffnessMatrices<TPZSkylineStructMatrix>(4);
}

BOOST_AUTO_TEST_SUITE_END()


TPZGeoMesh *threadTest::CreateGMesh(const int nDiv, int&matIdVol)
{
  constexpr MMeshType meshType = MMeshType::ETriangular;
  
  static TPZManVector<REAL,2> minX(2,0);
  static TPZManVector<REAL,2> maxX(2,1);
  TPZVec<int> nDivs(dim,nDiv);
  TPZManVector<int,1> matIdVec(1);
  matIdVol = matIdVec[0];
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(threadTest::dim,minX,maxX,matIdVec,nDivs,meshType,false);
}

TPZCompMesh *threadTest::CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol)
{
  auto *cmesh = new TPZCompMesh(gmesh);
  auto *laplacianMat = new TPZMatLaplacian(matIdVol, threadTest::dim);
  cmesh->InsertMaterialObject(laplacianMat);
 
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
  TPZGeoMesh *gMesh = CreateGMesh(nDiv, matIdVol);

  
  auto * cMesh = CreateCMesh(gMesh,pOrder,matIdVol);
    
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

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedSerial = end - start;
  std::cout << "Serial time: " << elapsedSerial.count() << "s\n";

  start = std::chrono::system_clock::now();
  auto matParallel = GetMatrix(nThreads);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedParallel = end - start;
  std::cout << "Parallel time: " << elapsedParallel.count() << "s\n";
  std::cout << "Speedup: " << elapsedSerial / elapsedParallel << "x\n";

  const int nr = matParallel->Rows();
  const int nc = matParallel->Cols();
  
  TPZFMatrix<STATE> matDiff(nr,nc, 0.0);
  matSerial->Substract(matParallel, matDiff);
  const auto normDiff = Norm(matDiff);
  std::cout.precision(17);
  std::cout << std::fixed << normDiff << std::endl;
  const bool checkMatNorm= IsZero(normDiff);
  BOOST_CHECK_MESSAGE(checkMatNorm,"failed");
  delete gMesh;
}
