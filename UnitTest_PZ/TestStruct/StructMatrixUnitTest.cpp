/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include "TPZGeoMeshTools.h"
#include "pzcmesh.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMaterialDataT.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzmatrix.h"
#include "TPZMatrixSolver.h"
#include "pzstepsolver.h"
#include "tpzautopointer.h"

//parallel layer classes
#include "pzstrmatrixor.h"
#include "pzstrmatrixot.h"
#ifdef PZ_USING_TBB
#include "pzstrmatrixflowtbb.h"
#endif
//struct matrices
#include "pzfstrmatrix.h"
#include "pzsfstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzsbstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"


#include <catch2/catch.hpp>
using namespace Catch::literals;


namespace structTest{
  /*!
    Creates computational mesh for dummy 1d problem
  */
  TPZAutoPointer<TPZCompMesh> CreateCMesh1D();
  TPZAutoPointer<TPZCompMesh> CreateCMesh2D();

  template <class TSTMAT>
  void CheckStiffnessMatrices(TPZAutoPointer<TPZCompMesh> cMesh,
                              const int nThreads);
  template <class TSTMAT>
  void CompareSerialParallelStiffMat(TPZAutoPointer<TPZCompMesh> cMesh,
                                     const int nThreads);
  template <class TSTMAT>
  void TestEquationFilter(TPZAutoPointer<TPZCompMesh> cMesh,
                                     const int nThreads);
  template <class TSTMAT1, class TSTMAT2>
  void CompareParallelLayerStiffMat(TPZAutoPointer<TPZCompMesh> cMesh,
                                     const int nThreads);
  
}

TEMPLATE_TEST_CASE("Assemble known matrix",
                   "[struct_tests][struct][multithread]",
                   TPZStructMatrixOR<STATE>,TPZStructMatrixOT<STATE>
#if defined(PZ_USING_TBB) && !defined(MACOSX)
                   ,TPZStructMatrixTBBFlow<STATE>
#endif
                   )
{
  auto cMesh = structTest::CreateCMesh1D();
  constexpr int nThreads{4};
  SECTION("Testing Full matrices"){
    SECTION("Non-symmetric"){
    structTest::CheckStiffnessMatrices<TPZFStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CheckStiffnessMatrices<TPZSFStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  SECTION("Testing Band matrices"){
    SECTION("Non-Symmetric"){
      structTest::CheckStiffnessMatrices<TPZBandStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CheckStiffnessMatrices<TPZSBandStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  SECTION("Testing Skyline matrices"){
    SECTION("Non-Symmetric"){
      structTest::CheckStiffnessMatrices<TPZSkylineNSymStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CheckStiffnessMatrices<TPZSkylineStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  SECTION("Testing Sparse matrices"){
    SECTION("Non-Symmetric"){
      structTest::CheckStiffnessMatrices<TPZSpStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CheckStiffnessMatrices<TPZSSpStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  cMesh = nullptr;
}



TEMPLATE_TEST_CASE("Compare parallel and serial matrices","[struct_tests][struct][multithread]",
                   TPZStructMatrixOR<STATE>,TPZStructMatrixOT<STATE>)
{
  auto cMesh = structTest::CreateCMesh2D();
  constexpr int nThreads{4};
  SECTION("Testing Full matrices"){
    SECTION("Non-symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZFStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZSFStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  SECTION("Testing Band matrices"){
    SECTION("Non-Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZBandStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZSBandStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  SECTION("Testing Skyline matrices"){
    SECTION("Non-Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZSkylineNSymStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZSkylineStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  SECTION("Testing Sparse matrices"){
    SECTION("Non-Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZSpStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
    SECTION("Symmetric"){
      structTest::CompareSerialParallelStiffMat<TPZSSpStructMatrix<STATE, TestType>>(cMesh, nThreads);
    }
  }
  cMesh = nullptr;
}

TEMPLATE_TEST_CASE("Test Equation Filter support",
                   "[struct_tests][struct][multithread]",
                   TPZStructMatrixOR<STATE>,TPZStructMatrixOT<STATE>
#if defined(PZ_USING_TBB) && !defined(MACOSX)
                   ,TPZStructMatrixTBBFlow<STATE>
#endif
                   )
{
  auto cMesh = structTest::CreateCMesh2D();
  constexpr int nThreads{4};

  structTest::TestEquationFilter<TestType>(cMesh, nThreads);
  cMesh = nullptr;
}

TEST_CASE("Compare parallel strategies","[struct_tests][struct][multithread]")
{
  auto cMesh = structTest::CreateCMesh2D();
  constexpr int nThreads{4};
  SECTION("Compare OR/OT"){
    structTest::CompareParallelLayerStiffMat<TPZStructMatrixOR<STATE>,TPZStructMatrixOT<STATE>>(cMesh, nThreads);
  }
#if defined(PZ_USING_TBB) && !defined(MACOSX)
  SECTION("Compare OR/TBBFlow"){
    structTest::CompareParallelLayerStiffMat<TPZStructMatrixOR<STATE>,TPZStructMatrixTBBFlow<STATE>>(cMesh, nThreads);
  }
  SECTION("Compare OT/TBBFlow"){
    structTest::CompareParallelLayerStiffMat<TPZStructMatrixOT<STATE>,TPZStructMatrixTBBFlow<STATE>>(cMesh, nThreads);
  }
#endif
  cMesh = nullptr;
}

namespace structTest{
  template <class TSTMAT>
  void CheckStiffnessMatrices(TPZAutoPointer<TPZCompMesh> cMesh,
                              const int nThreads)
  {

    TPZAutoPointer<TPZMatrix<STATE>> mat = [&cMesh,nThreads](){
      constexpr bool optimizeBandwidth{false};
      TPZLinearAnalysis an(cMesh, optimizeBandwidth);
      TSTMAT matskl(cMesh);
      matskl.SetNumThreads(nThreads);
      an.SetStructuralMatrix(matskl);
      TPZStepSolver<STATE> defaultSolver;
      defaultSolver.SetDirect(ELU);
      an.SetSolver(defaultSolver);
      an.Assemble();
      return an.MatrixSolver<STATE>().Matrix();
    }();
    const int nEq = cMesh->NElements() + 1;
    
    // mat->Print(std::cout);
    auto oldPrecision = Catch::StringMaker<STATE>::precision;
    REQUIRE(mat->GetVal(0,0) == 2.0_a);
    REQUIRE(mat->GetVal(nEq-1,nEq-1) == 2.0_a);
    REQUIRE(mat->GetVal(0,1) == 1.0_a);
    REQUIRE(mat->GetVal(nEq-1,nEq-2) == 1.0_a);
    for(auto i = 1; i < nEq-1; i++){
      REQUIRE(mat->GetVal(i,i-1) == 1.0_a);
      REQUIRE(mat->GetVal(i,i) == 4.0_a);
      REQUIRE(mat->GetVal(i,i+1) == 1.0_a);
    }
    Catch::StringMaker<STATE>::precision = oldPrecision;
  }

  template <class TSTMAT>
  void CompareSerialParallelStiffMat(TPZAutoPointer<TPZCompMesh> cMesh,
                                     const int nThreads)
  {

    auto GetMatrix = [cMesh](const int nThreads){
      constexpr bool optimizeBandwidth{false};
      TPZLinearAnalysis an(cMesh, optimizeBandwidth);
      TSTMAT matskl(cMesh);
      matskl.SetNumThreads(nThreads);
      an.SetStructuralMatrix(matskl);
      TPZStepSolver<STATE> defaultSolver;
      defaultSolver.SetDirect(ELU);
      an.SetSolver(defaultSolver);
      an.Assemble();
      return an.MatrixSolver<STATE>().Matrix();
    };
    auto matSerial = GetMatrix(0);
    auto matParallel = GetMatrix(nThreads);
    
    const int nr = matParallel->Rows();
    const int nc = matParallel->Cols();
  
    TPZFMatrix<STATE> matDiff(nr,nc, 0.0);
    matSerial->Subtract(matParallel, matDiff);
    const auto normDiff = Norm(matDiff);
    auto oldPrecision = Catch::StringMaker<STATE>::precision;
    CAPTURE(normDiff);
    REQUIRE(normDiff == Approx(0.0).margin(
                10*std::numeric_limits<STATE>::epsilon()));
    Catch::StringMaker<STATE>::precision = oldPrecision;
  }

  template <class TPARLAYER>
  void TestEquationFilter(TPZAutoPointer<TPZCompMesh> cMesh,
                                     const int nThreads)
  {

    auto GetMatrix = [cMesh](const int nThreads){
      constexpr bool optimizeBandwidth{false};
      TPZLinearAnalysis an(cMesh, optimizeBandwidth);
      const auto neqOld = cMesh->NEquations();
      const int neq = neqOld/2;
      TPZVec<int64_t> activeEqs(neq);
      for(int i = 0; i < neq; i++){
        activeEqs[i] = i;
      }

      TPZSkylineStructMatrix<STATE,TPARLAYER> matskl(cMesh);
      matskl.SetNumThreads(nThreads);
      matskl.EquationFilter().SetActiveEquations(activeEqs);
      an.SetStructuralMatrix(matskl);
      TPZStepSolver<STATE> defaultSolver;
      defaultSolver.SetDirect(ELU);
      an.SetSolver(defaultSolver);
      an.Assemble();
      return an.MatrixSolver<STATE>().Matrix();
    };
    auto matParallel = GetMatrix(nThreads);

    CAPTURE(true);

#ifdef PZ_USING_TBB    
    if constexpr (std::is_same_v<TPARLAYER, TPZStructMatrixTBBFlow<STATE>>){
      return;
    }
#endif
    auto matSerial = GetMatrix(0);
    CAPTURE(true);
  }


  template <class TSTMAT1, class TSTMAT2>
  void CompareParallelLayerStiffMat(TPZAutoPointer<TPZCompMesh> cMesh,
                                     const int nThreads)
  {
    
    auto GetMatrix = [cMesh](const int nThreads, bool first){
      constexpr bool optimizeBandwidth{false};
      TPZLinearAnalysis an(cMesh, optimizeBandwidth);
      TPZAutoPointer<TPZStructMatrix>strmat;
      if(first){
        strmat = new TPZSkylineStructMatrix<STATE, TSTMAT1>(cMesh);
      }
      else{
        strmat = new TPZSkylineStructMatrix<STATE,TSTMAT2>(cMesh);
      }
      strmat->SetNumThreads(nThreads);
      an.SetStructuralMatrix(strmat);
      TPZStepSolver<STATE> defaultSolver;
      defaultSolver.SetDirect(ELU);
      an.SetSolver(defaultSolver);
      an.Assemble();
      return an.MatrixSolver<STATE>().Matrix();
    };
    auto mat1 = GetMatrix(nThreads,true);
    auto mat2 = GetMatrix(nThreads,false);
    
    const int nr = mat2->Rows();
    const int nc = mat2->Cols();
  
    TPZFMatrix<STATE> matDiff(nr,nc, 0.0);
    mat1->Subtract(mat2, matDiff);
    const auto normDiff = Norm(matDiff);
    auto oldPrecision = Catch::StringMaker<STATE>::precision;
    CAPTURE(normDiff);
    REQUIRE(normDiff == Approx(0.0).margin(
                10*std::numeric_limits<STATE>::epsilon()));
    Catch::StringMaker<STATE>::precision = oldPrecision;
  }
  

  /****************************************
   *         AUXILIARY CLASSES            *
   *           AND FUNCTIONS              *
   ****************************************/

  
  class TPZMatTest : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>{
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>;
    const int fDim{0};
  public:
    TPZMatTest(const int matId, const int dim) :
      TBase(matId), fDim(dim) {}
                                                 
    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    int Dimension() const override{return 1;}
    int NStateVariables() const override{return 1;}
    void ContributeBC(const TPZMaterialDataT<STATE> &, REAL, TPZFMatrix<STATE>&,
                      TPZFMatrix<STATE>&,TPZBndCondT<STATE>&) override{;}
  };

  void TPZMatTest::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    if(fDim==1){
      ek.PutVal(0,0,2);
      ek.PutVal(1,1,2);
      ek.PutVal(0,1,1);
      ek.PutVal(1,0,1);
    }else{
      const TPZFMatrix<REAL>  &phi = data.phi;
      const auto phr = phi.Rows();
      for( auto in = 0; in < phr; in++ ) {
        for( auto jn = 0; jn < phr; jn++ ) {
          ek(in,jn) = phi(in,0) * phi(jn,0);
        }
      }
    }
  }

  TPZAutoPointer<TPZCompMesh> CreateCMesh1D()
  {
    constexpr int nEl{20};
    const TPZManVector<int,1> matIdVec={1};
    constexpr bool createBoundEls{false};
    TPZAutoPointer<TPZGeoMesh> gmesh =
      TPZGeoMeshTools::CreateGeoMesh1D(0, nEl, nEl,
                                       matIdVec, createBoundEls);
	constexpr int dim{1};
	TPZMatTest *material = new TPZMatTest(matIdVec[0],dim);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    constexpr int pOrder{1};
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(material);
	cmesh->AutoBuild();	
	return cmesh;
  }


  TPZAutoPointer<TPZCompMesh> CreateCMesh2D()
  {
    constexpr MMeshType meshType = MMeshType::ETriangular;
    constexpr int dim{2};
    constexpr int nDiv{4};
    constexpr bool createBoundEls{false};
    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,0};
    TPZManVector<int,2> nDivs = {nDiv,nDiv};
    TPZManVector<int,5> matIdVec = {1};
    TPZAutoPointer<TPZGeoMesh> gmesh =
      TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,
                                           nDivs,meshType,
                                           createBoundEls);
    
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    TPZMatTest *material = new TPZMatTest(matIdVec[0],dim);
    constexpr int pOrder{2};
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(material);
	cmesh->AutoBuild();	
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
  }

}
