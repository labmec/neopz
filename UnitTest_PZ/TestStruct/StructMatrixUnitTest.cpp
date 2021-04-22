/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"
#include "pzanalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
#include "tpzautopointer.h"

//parallel layer classes
#include "pzstrmatrixor.h"
#include "pzstrmatrixot.h"

//struct matrices
#include "pzfstrmatrix.h"
#include "pzsfstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzsbstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzbdstrmatrix.h"
#include "TPZBSpStructMatrix.h"


#include <catch2/catch.hpp>
using namespace Catch::literals;


namespace structChkTest{
  /*!
    Creates a 1D mesh with nel elements and size nel.
    Each element has unitary length
  */
  TPZAutoPointer<TPZGeoMesh>CreateGMesh1D(int64_t nel);
  /*!
    Creates computational mesh for dummy 1d problem
  */
  TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh>gmesh);

  template <class TSTMAT1>
  void CheckStiffnessMatrices(const int nThreads);
}


namespace structCprTest {
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
  constexpr int nThreads{4};
  SECTION("Testing Full matrices"){
    structChkTest::CheckStiffnessMatrices<TPZFStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZFStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);

    structChkTest::CheckStiffnessMatrices<TPZSFStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZSFStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);
  }
  SECTION("Testing Band matrices"){
    structChkTest::CheckStiffnessMatrices<TPZBandStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZBandStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);

    structChkTest::CheckStiffnessMatrices<TPZSBandStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZSBandStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);
  }
  SECTION("Testing Skyline matrices"){
    structChkTest::CheckStiffnessMatrices<TPZSkylineNSymStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZSkylineNSymStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);

    structChkTest::CheckStiffnessMatrices<TPZSkylineStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZSkylineStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);
  }
  SECTION("Testing Sparse matrices"){
    // structChkTest::CheckStiffnessMatrices<TPZSpStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    // structChkTest::CheckStiffnessMatrices<TPZSpStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);

    structChkTest::CheckStiffnessMatrices<TPZSSpStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    structChkTest::CheckStiffnessMatrices<TPZSSpStructMatrix<STATE, TPZStructMatrixOT<STATE>>>(nThreads);
  }
  // SECTION("Testing Block matrices"){
    // structChkTest::CheckStiffnessMatrices<TPZBlockDiagonalStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
    // structChkTest::CheckStiffnessMatrices<TPZBSpStructMatrix<STATE, TPZStructMatrixOR<STATE>>>(nThreads);
  // }
}

namespace structChkTest{
  /*!
    Creates a 1D mesh with nel elements and size nel.
    Each element has unitary length
  */
  TPZAutoPointer<TPZGeoMesh>CreateGMesh1D(int64_t nel)
  {
    TPZAutoPointer<TPZGeoMesh>gmesh = new TPZGeoMesh;
    constexpr int matId{1};
    const auto nnodes = nel + 1;
    gmesh->NodeVec().Resize(nnodes);
    for (auto i = 0 ; i < nnodes; i++)
      {
        const REAL pos = i;
        TPZManVector <REAL,3> coord= {pos,0.,0.};
        coord[0] = pos;
        gmesh->NodeVec()[i].SetCoord(coord);
        gmesh->NodeVec()[i].SetNodeId(i);
      }
    TPZManVector <int64_t,2> nodeVec(2);
    int64_t id;
    for (auto iel = 0; iel < nel; iel++)
      {
        nodeVec={iel,iel+1};
        gmesh->CreateGeoElement(EOned, nodeVec, matId, id);
        gmesh->ElementVec()[id];
      }
    gmesh->BuildConnectivity();
    return gmesh;
  }

  class TPZMatTest : public TPZMaterial{
  public:
    using TPZMaterial::TPZMaterial;
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    int Dimension() const override{return 1;}
    int NStateVariables() const override{return 1;}
    void ContributeBC(TPZMaterialData &, REAL, TPZFMatrix<STATE>&,
                      TPZFMatrix<STATE>&,TPZBndCond&) override{;}
  };

  void TPZMatTest::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    ek.PutVal(0,0,2);
    ek.PutVal(1,1,2);
    ek.PutVal(0,1,1);
    ek.PutVal(1,0,1);
  }

  TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh>gmesh)
  {
	constexpr int dim{1};
	constexpr int matId{1};
	TPZMatTest *material = new TPZMatTest(matId);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    constexpr int pOrder{1};
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(material);
	cmesh->AutoBuild();	
	return cmesh;
  }


  template <class TSTMAT>
  void CheckStiffnessMatrices(const int nThreads)
  {
    constexpr int nEl{20};
    auto gMesh = CreateGMesh1D(nEl);
    auto cMesh = CreateCMesh(gMesh);
    constexpr bool optimiseBandwidth{false};

    auto mat = [cMesh,nThreads](){
      constexpr bool optimizeBandwidth{false};
      TPZAnalysis an(cMesh, optimizeBandwidth);
      TSTMAT matskl(cMesh);
      matskl.SetNumThreads(nThreads);
      an.SetStructuralMatrix(matskl);
      an.Assemble();
      return an.MatrixSolver<STATE>().Matrix();
    }();

    constexpr int nEq{nEl+1};
    mat->Print(std::cout);
    REQUIRE(mat->GetVal(0,0) == 2.0_a);
    REQUIRE(mat->GetVal(nEq-1,nEq-1) == 2.0_a);
    REQUIRE(mat->GetVal(0,1) == 1.0_a);
    REQUIRE(mat->GetVal(nEq-1,nEq-2) == 1.0_a);
    for(auto i = 1; i < nEq-1; i++){
      REQUIRE(mat->GetVal(i,i-1) == 1.0_a);
      REQUIRE(mat->GetVal(i,i) == 4.0_a);
      REQUIRE(mat->GetVal(i,i+1) == 1.0_a);
    }
  }
}