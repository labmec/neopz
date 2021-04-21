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
#include "TPZMatLaplacian.h"
#include "TPZGeoMeshTools.h"
#include "pzbndcond.h"

#include "pzstrmatrixor.h"
#include "pzstrmatrixot.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include <catch2/catch.hpp>


namespace structTest {
  constexpr int dim{2};//aux variable
  //aux function for creating 2d gmesh on unit square
  TPZGeoMesh *CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC);
  //aux function for creating 2d cmesh with laplacian mat
  TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC);

  //test the stiffness matrices in serial and parallel computations
  template <class TSTMAT1, class TSTMAT2>
  void CompareStiffnessMatrices(const int nThreads);
}


TEST_CASE("structmatrix_assemble_test","[struct_tests]")
{
  SECTION("Testing Skyline matrices"){ 
  structTest::CompareStiffnessMatrices<
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOT<STATE>>
                           >(4);
  structTest::CompareStiffnessMatrices<
    TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixOT<STATE>>
                           >(4);
#ifdef PZ_USING_BOOST
  structTest::CompareStiffnessMatrices<
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSkylineStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>
                           >(4);
  structTest::CompareStiffnessMatrices<
    TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>
                           >(4);
#endif
  }
  SECTION("Testing Sparse matrices"){ 
  structTest::CompareStiffnessMatrices<
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSSpStructMatrix<STATE,TPZStructMatrixOT<STATE>>
                           >(4);
  structTest::CompareStiffnessMatrices<
    TPZSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSpStructMatrix<STATE,TPZStructMatrixOT<STATE>>
                           >(4);
#ifdef PZ_USING_BOOST
  structTest::CompareStiffnessMatrices<
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>
                           >(4);
  structTest::CompareStiffnessMatrices<
    TPZSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>,
    TPZSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>
                           >(4);
#endif
  }
}


TPZGeoMesh *structTest::CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC)
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
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(structTest::dim,minX,maxX,matIdVec,nDivs,meshType,true);
}

TPZCompMesh *structTest::CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC)
{
  auto *cmesh = new TPZCompMesh(gmesh);
  auto *laplacianMat = new TPZMatLaplacian(matIdVol, structTest::dim);

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

template <class TSTMAT1,class TSTMAT2>
void structTest::CompareStiffnessMatrices(const int nThreads)
{
  constexpr int nDiv{4};
  constexpr int pOrder{2};
  int matIdVol;
  int matIdBC;
  auto *gMesh = CreateGMesh(nDiv, matIdVol, matIdBC);
  auto *cMesh = CreateCMesh(gMesh, pOrder, matIdVol, matIdBC);
  
  //lambda for obtaining the FE matrix
  enum EWhich{EFirst,ESecond};
  auto GetMatrix = [cMesh,nThreads](EWhich which){
    constexpr bool optimizeBandwidth{false};
    TPZAnalysis an(cMesh, optimizeBandwidth);
    TPZAutoPointer<TPZStructMatrix> matskl =
      [cMesh,which]()->TPZStructMatrix*{
      
        TPZStructMatrix * mat_skl = nullptr;
        switch(which){
        case EFirst: mat_skl = new TSTMAT1(cMesh);break;
        case ESecond: mat_skl = new TSTMAT2(cMesh);break;
        }
      return mat_skl;
      }();
    
    matskl->SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    an.Assemble();
    return an.MatrixSolver<STATE>().Matrix();
  };

  auto start = std::chrono::system_clock::now();
  auto mat1 = GetMatrix(EFirst);
  std::cout<<typeid(TSTMAT1).name()<<std::endl;
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed1 = end - start;
  std::cout << "\tTime: " << elapsed1.count() << "s\n";

  start = std::chrono::system_clock::now();
  auto mat2 = GetMatrix(ESecond);
  std::cout<<typeid(TSTMAT2).name()<<std::endl;
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed2 = end - start;
  std::cout << "\tTime: " << elapsed2.count() << "s\n";

  const int nr = mat1->Rows();
  const int nc = mat1->Cols();
  
  TPZFMatrix<STATE> matDiff(nr,nc, 0.0);
  mat1->Substract(mat2, matDiff);
  const auto normDiff = Norm(matDiff);
  std::cout.precision(17);
  std::cout << std::fixed << "\tNorm diff: " << normDiff << std::endl;
  const bool checkMatNorm = IsZero(normDiff);
  REQUIRE(checkMatNorm);
  delete cMesh;
  delete gMesh;
}