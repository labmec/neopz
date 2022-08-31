//
//  TestSideShapeF.cpp
//  PZ
//
//  Created by Jeferson Fernandes on 31/08/22.
//  Copyright 2022 UNICAMP. All rights reserved.
//



#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapetetra.h"
#include "tpzquadrilateral.h"
#include "pzgeotetrahedra.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
using namespace pzshape;

#include "TPZShapeHDiv.h"
#include "pzelchdiv.h"

#include "TPZGeoMeshTools.h"
#include "TPZNullMaterial.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMaterialDataT.h"

#define USE_MAIN

#ifndef USE_MAIN
#include<catch2/catch.hpp>
#endif

using namespace pztopology;

enum SpaceType {EH1, EHDivConstant, EHDivKernel, EHDivStandard, EHCurl, EHCurlNoGrads};


TPZCompMesh * CreateCMeshH1(TPZGeoMesh* gmesh, int pOrder, const int volId);
TPZCompMesh * CreateCMeshHDiv(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam);

/*
    Test SideShapeFunction problem
*/
template<class top, class tshape>
void TesteSideShapeFunction(const int &pOrder, SpaceType stype);

// Test Hanging Nodes: FOR DEBUGGING PURPOSES
#ifndef USE_MAIN
TEST_CASE("Side Shape Function", "[side_shape_test]") {
    std::cout << "Testing Side Shape Functions \n";

    const int pOrder = GENERATE(1);
    SpaceType sType = GENERATE(EHDivStandard);
    // SpaceType sType = GENERATE(EHDivConstant);
    // SpaceType sType = GENERATE(EHCurl);
    
    TesteSideShapeFunction<pzshape::TPZShapeTriang>(pOrder,sType);
    // TesteSideShapeFunction<pzshape::TPZShapeQuad>(pOrder,sType);
    // TesteSideShapeFunction<pzshape::TPZShapeTetra>(pOrder,sType);
    // TesteSideShapeFunction<pzshape::TPZShapeCube>(pOrder,sType);
    std::cout << "Finish test Side Shape Functions \n";
}
#else
int main(){

    const int pOrder = 1;
    // SpaceType sType = EHCurl;
    SpaceType sType = EHDivStandard;
    // SpaceType sType = EHCurlNoGrads;
    // SpaceType sType = EHDivConstant;

    TesteSideShapeFunction<TPZQuadrilateral,pzshape::TPZShapeQuad>(pOrder,sType);

    return 0;
}
#endif

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template<class top, class tshape>
void TesteSideShapeFunction(const int &pOrder, SpaceType stype){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    const  int64_t pOrderIntRule = 5;
    const auto nSides = top::NSides;
    const auto nCorner = top::NCornerNodes;
    const auto dim = top::Dimension;
    const auto nFaces = top::NFacets;
    auto type = top::Type();
    const REAL tol = 1.e-6;
    const int volId = 1;
    const int BCId = 2;
    TPZGeoMesh *gmesh = TPZGeoMeshTools::CreateGeoMeshSingleElT<top>(volId,true,BCId);

#ifdef PZDEBUG
    //Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
#endif


    TPZCompMesh *cmesh;

    switch (stype)
    {
    case EH1:
        cmesh = CreateCMeshH1(gmesh, pOrder, volId);
        break;

    case EHDivStandard:
        cmesh = CreateCMeshHDiv(gmesh,pOrder,volId,HDivFamily::EHDivStandard);
        break;

    case EHDivConstant:
        cmesh = CreateCMeshHDiv(gmesh,pOrder,volId,HDivFamily::EHDivConstant);
        break;

    case EHDivKernel:
        cmesh = CreateCMeshHDiv(gmesh,pOrder,volId,HDivFamily::EHDivKernel);
        break;
            
    case EHCurl:
        // cmesh = CreateCMeshHCurl(gmesh,pOrder,volId,HCurlFamily::EHCurlStandard);
        // break;

    case EHCurlNoGrads:
        // cmesh = CreateCMeshHCurl(gmesh,pOrder,volId,HCurlFamily::EHCurlNoGrads);
        // break;

    default:
        break;
    };


    TPZCompElHDiv<tshape> *cel = dynamic_cast<TPZCompElHDiv<tshape> *>(cmesh->ElementVec()[0]);
    TPZMaterialDataT<STATE> data; 
    cel->InitMaterialData(data);
    
    TPZGeoEl *gel = cel->Reference();
    

    for (int iFace = nSides-nFaces-1; iFace < nSides-1; iFace++)
    {
        TPZIntPoints* intRule = gel->CreateSideIntegrationRule(iFace, pOrderIntRule);   
        TPZVec<REAL> node(dim);
        TPZVec<REAL> nodeSide(dim-1);
        const int npts = intRule->NPoints();

        TPZCompElSide thisside(cel, iFace);
        TPZGeoElSide thisgeoside = thisside.Reference();
        
        TPZCompElSide largecompside(cel, nSides-1);
        TPZGeoElSide largegeoside = largecompside.Reference();

        TPZTransform<> t = thisgeoside.SideToSideTransform(largegeoside);

        for (auto ipt = 0; ipt < npts; ipt++) 
        {
            REAL w;
            intRule->Point(ipt, nodeSide, w);
            
            t.Apply(nodeSide, node);
            std::cout << "Node = " << node << std::endl;
            std::cout << "Node2 = " << nodeSide << std::endl;

            cel->ComputeRequiredData(data,node);
            // int numshape = TPZShapeHDiv<tshape>::ComputeNConnectShapeF(iFace,pOrder);
            int numshapeSide = cel->NConnectShapeF(iFace-(nSides-nFaces-1),pOrder);
            TPZFNMatrix<100, REAL> phis(numshapeSide, 1), dphis(dim, numshapeSide);
            cel->SideShapeFunction(iFace,nodeSide,phis,dphis);

            // int numshapeVol = cel->NConnectShapeF(nSides-1,pOrder);
            // TPZFNMatrix<100, REAL> phiV(numshapeVol, dim), dphiV(dim, numshapeVol);
            
            cel->ComputeShape(node,data);
            // TPZShapeHDiv<tshape>::Shape(node, data, phiV, dphiV);
            std::cout << "PhiS = " << phis << std::endl;
            std::cout << "PhiV = " << data.fDeformedDirections << std::endl;
            std::cout << "end! " << std::endl;
        }    
    }
    
    
}




// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * CreateCMeshH1(TPZGeoMesh* gmesh, int pOrder, const int volId){

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(pOrder);

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(volId,gmesh->Dimension());
    cmesh->InsertMaterialObject(mat);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * CreateCMeshHDiv(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam){
    //Create flux mesh
    TPZCompMesh *cmeshflux = new TPZCompMesh(gmesh);
    cmeshflux->SetDimModel(gmesh->Dimension());
    cmeshflux->SetDefaultOrder(pOrder);
    
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(volId,gmesh->Dimension());
    cmeshflux->InsertMaterialObject(mat);

    cmeshflux->ApproxSpace().SetHDivFamily(hdivfam);
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(gmesh->Dimension());
    cmeshflux->AutoBuild();    

    return cmeshflux;
}