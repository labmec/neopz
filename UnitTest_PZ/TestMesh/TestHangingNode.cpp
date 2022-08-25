//
//  TestHangingNode.cpp
//  PZ
//
//  Created by Jeferson Fernandes and Nathan Shauer on 25/08/22.
//  Copyright 2022 UNICAMP. All rights reserved.
//



#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

using namespace pzshape;

#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "Projection/TPZL2ProjectionHDiv.h" //for BC in a single point
#include "Projection/TPZL2ProjectionCS.h"
#include <TPZGeoMeshTools.h>
#include "TPZEnumApproxFamily.h"
#include "TPZLinearAnalysis.h"
#include "TPZVTKGenerator.h"
#include "TPZNullMaterial.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzlog.h"

// #define USE_MAIN

#ifndef USE_MAIN
#include<catch2/catch.hpp>

#endif

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZGeoMesh* CreateGeoMesh(const TPZVec<int> &nDivs, const int volId, const int bcId);


enum SpaceType {EH1, EHDivConstant, EHDivKernel, EHDivStandard, EHCurl};

/*
    Test KernelHdiv problem
*/
template<class tshape>
void TestConstrainedSpace(const int &xdiv, const int &pOrder, SpaceType stype);

void Refinement(TPZGeoMesh *gmesh);

TPZCompMesh * CreateCMeshH1(TPZGeoMesh* gmesh, int pOrder, const int volId);
TPZCompMesh * CreateCMeshHDiv(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam);
TPZCompMesh * CreateCMeshHDiv2(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam);

// Test Hanging Nodes: FOR DEBUGGING PURPOSES
#ifndef USE_MAIN
TEST_CASE("Constrained Space", "[constrained_space_test]") {
  std::cout << "Testing Hanging Nodes \n";

  const int xdiv = GENERATE(2);
  const int pOrder = GENERATE(1);
  SpaceType sType = GENERATE(EHDivConstant);

//   TestConstrainedSpace<pzshape::TPZShapeTriang>(xdiv,pOrder,sType);
  TestConstrainedSpace<pzshape::TPZShapeQuad>(xdiv,pOrder,sType);
//   TestConstrainedSpace<pzshape::TPZShapeTetra>(xdiv,pOrder,sType);
//   TestConstrainedSpace<pzshape::TPZShapeCube>(xdiv, pOrder, sType);
  std::cout << "Finish test Constrained Space \n";
}
#else
int main(){

    const int xdiv = 2;
    const int pOrder = 1;
    SpaceType sType = EHDivStandard;

    TestConstrainedSpace<pzshape::TPZShapeQuad>(xdiv,pOrder,sType);

    return 0;
}
#endif

//Create 
template<class tshape>
TPZGeoMesh* CreateGeoMesh(const TPZVec<int> &nDivs, const int volId, const int bcId)
{
    MMeshType meshType;

    switch (tshape::Type())
    {
    case ETriangle:
        meshType = MMeshType::ETriangular;
        break;
    case EQuadrilateral:
        meshType = MMeshType::EQuadrilateral;
        break;
    case ETetraedro:
        meshType = MMeshType::ETetrahedral;
        break;
    case ECube:
        meshType = MMeshType::EHexahedral;
        break;
        case EPrisma:
        meshType = MMeshType::EPrismatic;
        break;
    default:
        DebugStop();
    }


    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*tshape::Dimension+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(tshape::Dimension, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);

    return gmesh;
    
}

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = 1.;
};

template<class tshape>
void TestConstrainedSpace(const int &xdiv, const int &pOrder, SpaceType stype){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    constexpr int volId{1};
    constexpr int bcId{2};
    int DIM = tshape::Dimension;

    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    auto gmesh = CreateGeoMesh<tshape>(nDivs, volId, bcId);

    // Refinement(gmesh);


    TPZCompMesh *cmesh;

    switch (stype)
    {
    case EH1:
        cmesh = CreateCMeshH1(gmesh, pOrder, volId);
        break;

    case EHDivStandard:
        cmesh = CreateCMeshHDiv2(gmesh,pOrder,volId,HDivFamily::EHDivStandard);
        break;

    case EHDivConstant:
        cmesh = CreateCMeshHDiv2(gmesh,pOrder,volId,HDivFamily::EHDivConstant);
        break;

    case EHDivKernel:
        cmesh = CreateCMeshHDiv2(gmesh,pOrder,volId,HDivFamily::EHDivKernel);
        break;

    default:
        break;
    };
    
    

    TPZLinearAnalysis an(cmesh);

    an.Assemble();
    an.Solve();
    

    TPZVec<std::string> fields = {"Solution"};

    std::set<int> matids;
    matids.insert(volId);
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences(); // compute integral in the multiphysics mesh
    TPZVec<STATE> vecint = cmesh->Integrate(fields[0], matids);

    std::cout << "Integral = " << vecint.size() <<  std::endl;
    for (int i = 0; i < vecint.size(); i++)
    {
        std::cout << vecint[i] << std::endl;   
    }
    
    const std::string plotfile = "solution";//sem o .vtk no final
    constexpr int vtkRes{0};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

    vtk.Do();

    REAL tol = 1.e-10;
#ifndef USE_MAIN
    REQUIRE(fabs(vecint[0]-1.)<tol);
#endif
    

}


void Refinement(TPZGeoMesh *gmesh){

    TPZManVector<TPZGeoEl*,10> children;
    gmesh->ElementVec()[0]->Divide(children);

    children[0]->Divide(children);
    

}

TPZCompMesh * CreateCMeshH1(TPZGeoMesh* gmesh, int pOrder, const int volId){

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(pOrder);

    TPZL2Projection<> *mat = new TPZL2Projection<>(volId,gmesh->Dimension());
    mat->SetForcingFunction(forcefunction,4);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}


TPZCompMesh * CreateCMeshHDiv(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam){
    //Create flux mesh
    TPZCompMesh *cmeshflux = new TPZCompMesh(gmesh);
    cmeshflux->SetDimModel(gmesh->Dimension());
    cmeshflux->SetDefaultOrder(pOrder);
    
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(volId);
    cmeshflux->InsertMaterialObject(mat);

    cmeshflux->ApproxSpace().SetHDivFamily(hdivfam);
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(gmesh->Dimension());
    cmeshflux->AutoBuild();

    //Create pressure cmesh
    TPZCompMesh *cmeshpres = new TPZCompMesh(gmesh);
    cmeshpres->SetDimModel(gmesh->Dimension());

    TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(volId);
    cmeshpres->InsertMaterialObject(mat2);

    switch (hdivfam)
    {
    case HDivFamily::EHDivStandard:
        cmeshpres->SetDefaultOrder(pOrder);
        cmeshpres->SetAllCreateFunctionsContinuous();
        cmeshpres->ApproxSpace().CreateDisconnectedElements(true);
        break;
    case HDivFamily::EHDivConstant:
        cmeshpres->SetDefaultOrder(0);
        cmeshpres->SetAllCreateFunctionsDiscontinuous();
        break;
    
    default:
        DebugStop();//HDivKernel nao tem malha de pressao
    }
    
    cmeshpres->AutoBuild();
    int ncon = cmeshpres->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmeshpres->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(1);
    }

    TPZManVector<TPZCompMesh *,2> meshvector = {cmeshflux,cmeshpres};

    //Create Multiphysics cmesh
    TPZMultiphysicsCompMesh* cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(gmesh->Dimension());

    auto *matL2 = new TPZL2ProjectionCS<>(volId,gmesh->Dimension());
    matL2->SetForcingFunction(forcefunction,4);
    cmesh->InsertMaterialObject(matL2);

    TPZManVector<int> active(meshvector.size(),1);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    

    return cmesh;

}

TPZCompMesh * CreateCMeshHDiv2(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam){
    //Create flux mesh
    TPZCompMesh *cmeshflux = new TPZCompMesh(gmesh);
    cmeshflux->SetDimModel(gmesh->Dimension());
    cmeshflux->SetDefaultOrder(pOrder);
    
    TPZL2ProjectionHDiv<> *mat = new TPZL2ProjectionHDiv<>(volId,gmesh->Dimension());
    mat->SetForcingFunction(forcefunction,4);
    cmeshflux->InsertMaterialObject(mat);

    cmeshflux->ApproxSpace().SetHDivFamily(hdivfam);
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(gmesh->Dimension());
    cmeshflux->AutoBuild();    

    return cmeshflux;

}