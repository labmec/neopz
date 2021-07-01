/**
 * @file
 * @brief Contains Unit Tests for checking compatibility of approx spaces
 */

#include <catch2/catch.hpp>
#include "TPZMatDeRhamH1.h"
#include "TPZMatDeRhamHCurl.h"
#include "TPZMatDeRhamHDiv.h"
#include "TPZMatDeRhamL2.h"
#include "TPZMatDeRhamH1HCurl.h"
#include "TPZNullMaterial.h"
#include "TPZGeoMeshTools.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZLinearAnalysis.h"
#include "pzfstrmatrix.h"
#include "TPZMatrixSolver.h"
#include "TPZElementMatrixT.h"
#include "MMeshType.h"

enum class ESpace{H1,HCurl,HDiv,L2};

static const std::map<ESpace, const char *> names({{ESpace::H1, "H1"},
                                                   {ESpace::HCurl, "HCurl"},
                                                   {ESpace::HDiv, "HDiv"},
                                                   {ESpace::L2, "L2"}});


/** @brief Compares the dimensions of the span of the operator associated
    with the left space against the kernel of the right space*/
template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckCompatibilityUniformMesh(int k);
/** @brief Checks if the functions of the right space span the range of the
    operator of the left space*/
template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckExactSequence(int kRight);

TEMPLATE_TEST_CASE("Dimension Compatibility", "[derham_tests]",
                   (typename std::integral_constant<int,2>),
                   (typename std::integral_constant<int,3>))
{
  //SVD requires LAPACK
#ifndef PZ_USING_LAPACK
  return;
#endif
  constexpr int dim = TestType::value;

  ESpace leftSpace = GENERATE(ESpace::H1,
                              ESpace::HCurl,
                              ESpace::HDiv);
  int k = GENERATE(1,2,3,4,5);
  SECTION(names.at(leftSpace)) {
    switch (leftSpace) {
    case ESpace::H1:
      CheckCompatibilityUniformMesh<ESpace::H1,ESpace::HCurl,dim>(k);
      //TODOFIX
      // if constexpr (dim == 2){
      //   CheckCompatibilityUniformMesh<ESpace::H1, ESpace::HDiv, dim>(k);
      // }
      break;
    case ESpace::HCurl:
      if constexpr (dim == 2){
        CheckCompatibilityUniformMesh<ESpace::HCurl, ESpace::L2, dim>(k);
      }
      //TODOFIX
      // else{
      //   CheckCompatibilityUniformMesh<ESpace::HCurl, ESpace::HDiv, dim>(k);
      // }
      break;
    case ESpace::HDiv:
      CheckCompatibilityUniformMesh<ESpace::HDiv,ESpace::L2,dim>(k);
      break;
    case ESpace::L2:
      break;
    }
  }
}

TEMPLATE_TEST_CASE("Inclusion", "[derham_tests]",
                   (typename std::integral_constant<int,2>),
                   (typename std::integral_constant<int,3>))
{
  //SVD requires LAPACK
#ifndef PZ_USING_LAPACK
  return;
#endif
  constexpr int dim = TestType::value;

  ESpace leftSpace = GENERATE(ESpace::H1,
                              ESpace::HCurl,
                              ESpace::HDiv);
  int k = GENERATE(1,2,3,4,5);
  SECTION(names.at(leftSpace)) {
    switch (leftSpace) {
    case ESpace::H1:
      CheckExactSequence<ESpace::H1,ESpace::HCurl,dim>(k);
      //TODOFIX
      // if constexpr (dim == 2){
      //   CheckExactSequence<ESpace::H1, ESpace::HDiv, dim>(k);
      // }
      break;
      //TODOFIX
    default:
      break;
    }
  }
}

/*******************************************************************************
 ***********************************AUX FUNCS***********************************
*******************************************************************************/
TPZAutoPointer<TPZGeoMesh> CreateGMesh(const int dim, const MMeshType elType,
                                       const int matId);

TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh> gMesh,
                                        TPZMaterial *mat, const int matId,
                                        const int k, ESpace space);

int CalcRank(const TPZFMatrix<STATE> & S, const STATE tol);

/*******************************************************************************
***********************************TEST FUNCS***********************************
*******************************************************************************/

template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckCompatibilityUniformMesh(int kRight) {
  const auto nameLeft = names.at(leftSpace);
  const auto nameRight = names.at(rightSpace);
  constexpr int matId = 1;
  auto elType = GENERATE(MMeshType::ETriangular,
                         MMeshType::EQuadrilateral,
                         MMeshType::ETetrahedral,
                         MMeshType::EHexahedral,
                         MMeshType::EPrismatic);
  const auto elName = MMeshType_Name(elType);
  const auto elDim = MMeshType_Dimension(elType);

  //TODOFIX
  if constexpr(leftSpace==ESpace::HCurl||
               rightSpace==ESpace::HCurl){
    if(elType!=MMeshType::ETriangular &&
       elType!=MMeshType::ETetrahedral){
      return;
    }
  }
  if(elDim == dim){

    int kLeft;
    auto *matLeft = [&kLeft,kRight]()
      -> TPZMaterial* {
      if constexpr (leftSpace == ESpace::H1) {
        kLeft = kRight + 1;
        return new TPZMatDeRhamH1(matId, dim);
      } else if constexpr (leftSpace == ESpace::HCurl) {
        kLeft = kRight + 1;
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

    auto gmesh = CreateGMesh(dim,elType,matId);
    auto cmeshL = CreateCMesh(gmesh, matLeft, matId, kLeft, leftSpace);
    auto cmeshR = CreateCMesh(gmesh, matRight, matId, kRight, rightSpace);

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
    static constexpr auto tol = std::numeric_limits<STATE>::epsilon()*100000;
    const int rankL = CalcRank(SL,tol);
    const int rankR = CalcRank(SR,tol);
    const int kerR = dimR-rankR;
      
    CAPTURE(kLeft,kRight,dimL,rankL,dimR,kerR,rankR);

    if constexpr (rightSpace==ESpace::L2)//for L2 the operator is -> 0
      REQUIRE(rankL == rankR);
    else
      REQUIRE(rankL == kerR);
    
    delete matLeft;
    delete matRight;
  }
}

template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckExactSequence(int kRight) {
  const auto nameLeft = names.at(leftSpace);
  const auto nameRight = names.at(rightSpace);
  constexpr int matId = 1;
  auto elType = GENERATE(MMeshType::ETriangular, MMeshType::EQuadrilateral,
                         MMeshType::ETetrahedral, MMeshType::EHexahedral,
                         MMeshType::EPrismatic);
  const auto elName = MMeshType_Name(elType);
  const auto elDim = MMeshType_Dimension(elType);

  // TODOFIX
  if constexpr (leftSpace == ESpace::HCurl || rightSpace == ESpace::HCurl) {
    if (elType != MMeshType::ETriangular && elType != MMeshType::ETetrahedral) {
      return;
    }
  }
  if (elDim == dim) {

    const int kLeft = [kRight]() {
      if constexpr (leftSpace == ESpace::H1) {
        return kRight + 1;
      } else if constexpr (leftSpace == ESpace::HCurl) {
        return kRight + 1;
      } else if constexpr (leftSpace == ESpace::HDiv) {
        return kRight;
      } else if constexpr (leftSpace == ESpace::L2) {
        return -1;
      }
    }();
    auto *matLeft = new TPZNullMaterial<STATE>(matId,dim,1);
    auto *matRight = new TPZNullMaterial<STATE>(matId,dim,1);

    CAPTURE(nameLeft, nameRight, elName);
    constexpr bool createBoundEls{false};
    TPZAutoPointer<TPZGeoMesh> gmesh =
        TPZGeoMeshTools::CreateGeoMeshSingleEl(elType, matId, createBoundEls);
    
    auto cmeshL = CreateCMesh(gmesh, matLeft, matId, kLeft, leftSpace);
    auto cmeshR = CreateCMesh(gmesh, matRight, matId, kRight, rightSpace);

    TPZAutoPointer<TPZCompMesh> cmeshMF = new TPZCompMesh(gmesh);

    auto *matMF = []() -> TPZMaterial * {
      if constexpr (rightSpace == ESpace::H1) {
        return nullptr;
      } else if constexpr (rightSpace == ESpace::HCurl) {
        return new TPZMatDeRhamH1HCurl(matId, dim);
      } else if constexpr (rightSpace == ESpace::HDiv) {
        return nullptr;
      } else if constexpr (rightSpace == ESpace::L2) {
        return nullptr;
      }
    }();
    if (!matMF)
      return;

    cmeshMF->InsertMaterialObject(matMF);
    cmeshMF->SetDimModel(dim);

    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild();
    cmeshMF->CleanUpUnconnectedNodes();

    TPZManVector<TPZCompMesh *, 2> meshVecIn = {cmeshL.operator->(),
                                                cmeshR.operator->()};

    TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn,
                                                 cmeshMF.operator->());

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();
    constexpr bool optimizeBandwidth{false};
    TPZLinearAnalysis analysis(cmeshMF, optimizeBandwidth);
    TPZFStructMatrix<STATE> strmtrx(cmeshMF);
    analysis.SetStructuralMatrix(strmtrx);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    
    analysis.Assemble();
    //it is a full matrix
    auto mfMatrix =
      dynamic_cast<TPZFMatrix<STATE>&> (
          analysis.MatrixSolver<STATE>().Matrix().operator*());

    const int rLeftEqs = cmeshR->NEquations();
    const int lLeftEqs = cmeshL->NEquations();
    TPZFMatrix<STATE> rightMatrix(rLeftEqs,rLeftEqs);
    mfMatrix.GetSub(lLeftEqs, lLeftEqs, rLeftEqs, rLeftEqs, rightMatrix);


    TPZFMatrix<STATE> Sfull, Sright;
    
    {
      TPZFMatrix<STATE> Udummy, VTdummy;
      mfMatrix.SVD(Udummy, Sfull, VTdummy, 'N', 'N');
      rightMatrix.SVD(Udummy, Sright, VTdummy, 'N', 'N');
    }
    
    static constexpr auto tol = std::numeric_limits<STATE>::epsilon()*100;
    auto fullDim = Sfull.Rows();
    auto fullRank = CalcRank(Sfull,tol);
    auto rightDim = Sright.Rows();
    auto rightRank = CalcRank(Sright,tol);


    // if(elType == MMeshType::ETriangular){
    //   const int postProcessResolution = 3;
    //   const std::string executionInfo = [&]() {
    //     std::string name("");
    //     name.append(elName);
    //     name.append(std::to_string(kRight));
    //     return name;
    //   }();

    //   const std::string plotfile =
    //       "solution" + executionInfo + ".vtk"; // where to print the vtk files
    //   TPZStack<std::string> scalnames, vecnames;
    //   if constexpr (leftSpace == ESpace::HDiv || leftSpace == ESpace::HCurl) {
    //     vecnames.Push("SolutionLeft"); // print the state variable
    //   } else {
    //     scalnames.Push("SolutionLeft"); // print the state variable
    //   }

    //   if constexpr (rightSpace == ESpace::HDiv || rightSpace == ESpace::HCurl) {
    //     vecnames.Push("SolutionRight"); // print the state variable
    //   } else {
    //     scalnames.Push("SolutionRight"); // print the state variable
    //   }

    //   analysis.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    //   TPZFMatrix<STATE> sol(cmeshMF->NEquations(), 1);
    //   sol.Zero();
    //   for (int i = 0; i < sol.Rows(); i++) {
    //     sol(i - 1 < 0 ? 0 : i - 1, 0) = 0;
    //     sol(i, 0) = 1;
    //     analysis.LoadSolution(sol);
    //     TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(
    //         meshVecIn, cmeshMF.operator->());
    //     analysis.PostProcess(postProcessResolution);
    //   }
    // }

    CAPTURE(kLeft,kRight,fullDim, fullRank,rightDim, rightRank);
    REQUIRE(fullRank == rightRank);
  }
}


/*******************************************************************************
 ***********************************AUX FUNCS***********************************
*******************************************************************************/
TPZAutoPointer<TPZGeoMesh> CreateGMesh(const int dim, const MMeshType elType,
                                       const int matId)
{
  constexpr int ndiv{2};
  TPZManVector<int, 3> nDivVec(dim, ndiv);
  TPZManVector<REAL, 3> minX({0, 0, 0});
  TPZManVector<REAL, 3> maxX({1, 1, 1});
  if (dim == 2)
    maxX[2] = 0.0;
  constexpr bool createBoundEls{false};
  TPZAutoPointer<TPZGeoMesh> gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(
      dim, minX, maxX, {matId}, nDivVec, elType, createBoundEls);
  return gmesh;
}

TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                                        TPZMaterial *mat, const int matId,
                                        const int k, ESpace space)
{

  TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
  const int nel = cmesh->NElements();
  cmesh->SetDefaultOrder(k);
  cmesh->SetDimModel(gmesh->Dimension());
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
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    break;
  }
  cmesh->AutoBuild();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}



int CalcRank(const TPZFMatrix<STATE> & S, const STATE tol){
  int rank = 0;
  const int dimMat = S.Rows();
  for (int i = 0; i < dimMat; i++) {
    rank += S.GetVal(i, 0) > tol ? 1 : 0;
  }
  return rank;
};