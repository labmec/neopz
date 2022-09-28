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

// #define USE_MAIN

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
template<class tshape>
void TesteSideShapeFunction(const int &pOrder, SpaceType stype);

// Test Hanging Nodes: FOR DEBUGGING PURPOSES
#ifndef USE_MAIN
TEST_CASE("Side Shape Function", "[side_shape_test]") {
    std::cout << "Testing Side Shape Functions \n";

    const int pOrder = GENERATE(1,2,3);
    SpaceType sType = GENERATE(EHDivConstant,EHDivStandard);
    
    TesteSideShapeFunction<pzshape::TPZShapeTriang>(pOrder,sType);
    TesteSideShapeFunction<pzshape::TPZShapeQuad>(pOrder,sType);
    TesteSideShapeFunction<pzshape::TPZShapeTetra>(pOrder,sType);
    TesteSideShapeFunction<pzshape::TPZShapeCube>(pOrder,sType);
    std::cout << "Finish test Side Shape Functions \n";
}
#else
int main(){

    const int pOrder = 1;
    // SpaceType sType = EHCurl;
    // SpaceType sType = EHDivStandard;
    // SpaceType sType = EHCurlNoGrads;
    SpaceType sType = EHDivConstant;
    // SpaceType sType = EH1;

    TesteSideShapeFunction<pzshape::TPZShapeQuad>(pOrder,sType);
    TesteSideShapeFunction<pzshape::TPZShapeTriang>(pOrder,sType);
    TesteSideShapeFunction<pzshape::TPZShapeTetra>(pOrder,sType);
    TesteSideShapeFunction<pzshape::TPZShapeCube>(pOrder,sType);

    return 0;
}
#endif

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template<class tshape>
void TesteSideShapeFunction(const int &pOrder, SpaceType stype){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    const  int64_t pOrderIntRule = 5;
    const auto nSides = tshape::NSides;
    const auto nCorner = tshape::NCornerNodes;
    const auto dim = tshape::Dimension;
    const auto nFaces = tshape::NFacets;
    auto type = tshape::Type();
    MMeshType meshtype;
    switch (type){
    case ETriangle: 
        meshtype = MMeshType::ETriangular;
        break;
    case EQuadrilateral: 
        meshtype = MMeshType::EQuadrilateral;
        break;
    case ETetraedro: 
        meshtype = MMeshType::ETetrahedral;
        break;
    case ECube: 
        meshtype = MMeshType::EHexahedral;
        break;
    
    default:
        break;
    }
    const REAL tol = 1.e-6;
    const int volId = 1;
    const int BCId = 2;
    TPZGeoMesh *gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshtype,volId,true,BCId);

#ifdef PZDEBUG
    //Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
#endif


    TPZCompMesh *cmesh;
    TPZIntelGen<tshape> *cel;

    switch (stype)
    {
    case EH1:
        cmesh = CreateCMeshH1(gmesh, pOrder, volId);
        cel = dynamic_cast<TPZIntelGen<tshape> *>(cmesh->ElementVec()[0]);
        break;

    case EHDivStandard:
        cmesh = CreateCMeshHDiv(gmesh,pOrder,volId,HDivFamily::EHDivStandard);
        cel = dynamic_cast<TPZCompElHDiv<tshape> *>(cmesh->ElementVec()[0]);
        break;

    case EHDivConstant:
        cmesh = CreateCMeshHDiv(gmesh,pOrder,volId,HDivFamily::EHDivConstant);
        cel = dynamic_cast<TPZCompElHDiv<tshape> *>(cmesh->ElementVec()[0]);
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
        DebugStop();
        break;
    };

    

    
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
            // std::cout << "Node = " << node << std::endl;
            // std::cout << "Node2 = " << nodeSide << std::endl;

            cel->ComputeRequiredData(data,node);

            int numshapeSide = cel->NConnectShapeF(iFace-(nSides-nFaces-1),pOrder);
            TPZFNMatrix<100, REAL> phis(numshapeSide, 1), dphis(dim, numshapeSide);
            cel->SideShapeFunction(iFace,nodeSide,phis,dphis);
            cel->ComputeShape(node,data);

            TPZFNMatrix<100, REAL> phiV(numshapeSide, dim);
            if (stype == EHDivConstant || stype == EHDivStandard){
                for (int k = 0; k < numshapeSide; k++)
                {
                    for (int d = 0; d < dim; d++)
                    {
                        phiV(k,d)=data.fDeformedDirections(d,numshapeSide*(iFace-nSides+nFaces+1)+k);
                    }
                }
            } else {
                std::cout << "Test not implemented for this space \n";
                DebugStop();
            }
            

            TPZVec<REAL> normalSide;
            thisgeoside.Normal(nodeSide,normalSide);

            TPZManVector<REAL> res(numshapeSide,0.);
            for (int i = 0; i < phiV.Rows(); i++){
                for (int j = 0; j < phiV.Cols(); j++){
                    res[i] += phiV(i,j)*normalSide[j];
                }   
            }
            
            if (dim == 3 && stype == EHDivConstant){
                bool orient = false;
                switch (type)
                {
                case ETetraedro:
                    if (iFace == 10 || iFace == 13) orient = true;
                    break;
                case ECube:
                    if (iFace == 20 || iFace == 23 || iFace == 24) orient = true;
                    break;
                
                default:
                    break;
                }
                if (orient){
                    for (int i = 1; i < phiV.Rows(); i++) res[i] *= -1.;
                }
                
            }
            // std::cout << "PhiS = " << phis << std::endl;
            // std::cout << "PhiV = " << data.fDeformedDirections << std::endl;
            // std::cout << "PhiV = " << phiV << std::endl;
            // std::cout << "res = " << res << std::endl;
        REAL tol = 1.e-8;
#ifndef USE_MAIN
        for (int i = 0; i < numshapeSide; i++)
        {
            REQUIRE(res[i]-phis(i,0) < tol);
        }      
#else
        for (int i = 0; i < numshapeSide; i++)
        {
            if(fabs(res[i]-phis(i,0)) > tol){
                std::cout << "Problem!!\n";
                std::cout << "Res = " << res << std::endl;
                std::cout << "phis = " << phis << std::endl;
                std::cout << "diff = " << res[i]-phis(i,0)<< std::endl;
                DebugStop(); 
            }
        }      
#endif
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