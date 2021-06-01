/**
 * @file
 * @brief Contains Unit Tests for checking compatibility of approx spaces
 */

#include <catch2/catch.hpp>
#include "TPZMatDeRhamH1.h"
#include "TPZMatDeRhamHCurl.h"
#include "TPZMatDeRhamHDiv.h"
#include "TPZMatDeRhamL2.h"
#include "TPZGeoMeshTools.h"
#include "pzcmesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZMatrixSolver.h"
#include "TPZElementMatrixT.h"
#include "MMeshType.h"

enum class ESpace{H1,HCurl,HDiv,L2};

static const std::map<ESpace, const char *> names({{ESpace::H1, "H1"},
                                                   {ESpace::HCurl, "HCurl"},
                                                   {ESpace::HDiv, "HDiv"},
                                                   {ESpace::L2, "L2"}});

template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckCompatibilityUniformMesh(int k);

TEMPLATE_TEST_CASE("Compatibility Uniform Mesh", "[derham_tests]",
                   (typename std::integral_constant<int,2>),
                   (typename std::integral_constant<int,3>))
{
  //SVD requires LAPACK
#ifndef PZ_USING_LAPACK
  return;
#endif
  return;
  constexpr int dim = TestType::value;
  ESpace leftSpace = GENERATE(ESpace::H1, ESpace::HCurl, ESpace::HDiv);
  int k = GENERATE(1);
  SECTION(names.at(leftSpace)) {
    switch (leftSpace) {
    case ESpace::H1:
      SECTION(names.at(ESpace::HCurl)) {
        CheckCompatibilityUniformMesh<ESpace::H1,ESpace::HCurl,dim>(k);
      }
      if constexpr (dim == 2){
        SECTION(names.at(ESpace::HDiv)) {
          CheckCompatibilityUniformMesh<ESpace::H1, ESpace::HDiv, dim>(k);
        }
      }
      break;
    case ESpace::HCurl:
      if constexpr (dim == 2){
        SECTION(names.at(ESpace::L2)) {
          CheckCompatibilityUniformMesh<ESpace::HCurl, ESpace::L2, dim>(k);
        }
      }else{
        SECTION(names.at(ESpace::HDiv)) {
          CheckCompatibilityUniformMesh<ESpace::HCurl, ESpace::HDiv, dim>(k);
        }
      }
      break;
    case ESpace::HDiv:
      SECTION(names.at(ESpace::L2)) {
        CheckCompatibilityUniformMesh<ESpace::HDiv,ESpace::L2,dim>(k);
      }
      break;
    case ESpace::L2:
      break;
    }
  }
}

template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckCompatibilityUniformMesh(int kRight) {
  const auto nameLeft = names.at(leftSpace);
  const auto nameRight = names.at(rightSpace);
  constexpr int matId = 1;
  auto elType = GENERATE(MMeshType::ETriangular,
                         MMeshType::EQuadrilateral,
                         MMeshType::ETetrahedral,
                         MMeshType::EPrismatic);
  const auto elName = MMeshType_Name(elType);
  const auto elDim = MMeshType_Dimension(elType);
  
  if(elDim != dim) return;
  
  SECTION(elName){

    int kLeft;
    auto *matLeft = [&kLeft,kRight]()
      -> TPZMaterial* {
      if constexpr (leftSpace == ESpace::H1) {
        kLeft = kRight + 1;
        return new TPZMatDeRhamH1(matId, dim);
      } else if constexpr (leftSpace == ESpace::HCurl) {
        kLeft = kRight;
        return new TPZMatDeRhamHCurl(matId,dim);
      } else if constexpr (leftSpace == ESpace::HDiv) {
        kLeft = kRight;
        return new TPZMatDeRhamHDiv(matId,dim);
      } else if constexpr (leftSpace == ESpace::L2) {
        return nullptr;
      }
    }();
    
    auto *matRight = []()
      -> TPZMaterial * {
      if constexpr (rightSpace == ESpace::H1) {
        return nullptr;
      } else if constexpr (rightSpace == ESpace::HCurl) {
        return new TPZMatDeRhamHCurl(matId,dim);
      } else if constexpr (rightSpace == ESpace::HDiv) {
        return new TPZMatDeRhamHDiv(matId,dim);
      } else if constexpr (rightSpace == ESpace::L2) {
        return new TPZMatDeRhamL2(matId,dim);
      }
    }();

    CAPTURE(nameLeft,nameRight,elName);
    REQUIRE((matLeft != nullptr && matRight != nullptr));
    
    constexpr int ndiv{1};
    TPZManVector<int, 3> nDivVec(dim, ndiv);
    TPZManVector<REAL, 3> minX({0, 0, 0});
    TPZManVector<REAL, 3> maxX({1, 1, 1});
    if constexpr (dim == 2)
      maxX[2] = 0.0;
    constexpr bool createBoundEls{false};
    TPZAutoPointer<TPZGeoMesh> gmesh =
            TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX,maxX,
                                                 {matId},nDivVec,elType,
                                                 createBoundEls);

    auto CreateCMesh = [](TPZAutoPointer<TPZGeoMesh> gmesh,
                          TPZMaterial *mat,
                          const int k, ESpace space){
      TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
      const int nel = cmesh->NElements();
      cmesh->SetDefaultOrder(k);
      cmesh->SetDimModel(dim);
      cmesh->InsertMaterialObject(mat);
      switch (space) {
      case ESpace::H1:
        cmesh->SetAllCreateFunctionsContinuous();
        break;
      case ESpace::HCurl:
        cmesh->SetAllCreateFunctionsHCurl();
        break;
      case ESpace::HDiv:
        cmesh->SetAllCreateFunctionsHDiv();
        break;
      case ESpace::L2:
        cmesh->SetAllCreateFunctionsDiscontinuous();
        break;
      }
      cmesh->AutoBuild();
      cmesh->CleanUpUnconnectedNodes();
      return cmesh;
    };

    
    auto cmeshL = CreateCMesh(gmesh,matLeft,kLeft,leftSpace);
    auto cmeshR = CreateCMesh(gmesh,matRight,kRight,rightSpace);

    TPZAutoPointer<TPZMatrix<STATE>> matPtrL = nullptr;
    TPZAutoPointer<TPZMatrix<STATE>> matPtrR = nullptr;

    //for debugging
    constexpr bool singleElement{true};

    if constexpr (singleElement){
      TPZCompEl* celL = cmeshL->Element(0);
      TPZCompEl* celR = cmeshR->Element(0);
      
      const int neqL = celL->NEquations();
      const int neqR = celR->NEquations();
      
      TPZElementMatrixT<STATE> matKL, matKR, matFL, matFR;

      celL->CalcStiff(matKL, matFL);

      celR->CalcStiff(matKR, matFR);
      
      auto *matL = new TPZFMatrix<STATE>(matKL.Matrix());
      matPtrL = TPZAutoPointer<TPZMatrix<STATE>>(matL);

      auto *matR = new TPZFMatrix<STATE>(matKR.Matrix());
      matPtrR = TPZAutoPointer<TPZMatrix<STATE>>(matR);
    }else{
      TPZLinearAnalysis anL(cmeshL,false);
      TPZLinearAnalysis anR(cmeshR, false);
      anL.Assemble();
      anR.Assemble();

      matPtrL = anL.MatrixSolver<STATE>().Matrix();
      matPtrR = anR.MatrixSolver<STATE>().Matrix();
    }

    TPZFMatrix<STATE> matL(*matPtrL);
    TPZFMatrix<STATE> matR(*matPtrR);

    const int dimL = matL.Rows();
    const int dimR = matR.Rows();
    
    TPZFMatrix<STATE> SL, SR;
    
    {
      TPZFMatrix<STATE> Udummy, VTdummy;
      matL.SVD(Udummy, SL, VTdummy, 'N', 'N');
      matR.SVD(Udummy, SR, VTdummy, 'N', 'N');
    }

    constexpr auto tol = std::numeric_limits<STATE>::epsilon();
    auto CalcRank = [tol](const TPZFMatrix<STATE> & S){
      int rank = 0;
      const int dimMat = S.Rows();
      for(int i = 0; i < dimMat; i++){
        rank += S.GetVal(i,0) > tol ? 1 : 0;
      }
      return rank;
    };
    const int rankL = CalcRank(SL);
    const int rankR = CalcRank(SR);
    const int kerR = dimR-rankR;

    if(
        !strcmp(nameLeft,"HCurl")||
        !strcmp(nameRight,"HCurl")){
      ;//for now let us debug without hcurl
    }
    else{
      CAPTURE(kLeft,kRight,dimR,dimL,rankL,kerR,rankR);
      REQUIRE(rankL >= kerR);
    }
    delete matLeft;
    delete matRight;
  }
}