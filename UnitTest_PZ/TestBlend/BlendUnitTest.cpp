/**
 * @file BlendUnitTest.cpp
 * @brief Define a Unit Test using Boost for blend elements over all kind of geometries available in the library
 *
 */
#include <iostream>

#include "TPZRefPatternDataBase.h"

#include "pzgeotriangle.h"
#include "pzgeoquad.h"

#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"

#include "pzgeoelrefless.h"

#include "tpzquadraticline.h"

#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"

#include "tpzquadratictetra.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticpyramid.h"

#include "tpzgeoblend.h"

#include "TPZVTKGeoMesh.h"

#include "tpzgeoelrefpattern.h"

#include "tpzarc3d.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "TPZGenSpecialGrid.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testgeom");
#endif

#include "fad.h"


#include<catch2/catch.hpp>

//#define BLEND_VERBOSE //outputs x and grad comparisons
#define BLEND_OUTPUT_TXT//prints all elements in .txt format
#define BLEND_OUTPUT_VTK//prints all elements in .vtk format

namespace blendtest{
const int pOrder = 10;
const REAL tol = 1e-8;
//check the norm of the difference between two given vectors
bool CheckVectors(const TPZVec<REAL> &x1, std::string name1,
                  const TPZVec<REAL> &x2, std::string name2, const REAL &tol);
//check the fobrenius norm of the difference between two given matrices
bool CheckMatrices(const TPZFMatrix<REAL> &gradx1, std::string name1,
                   const TPZFMatrix<REAL> &gradx2, std::string name2, const REAL &tol);
template <class TGeo>
void CompareQuadraticAndBlendEls();
void TestMesh2D(TPZGeoMesh *gmesh, int nDiv);
void TestMesh3D(TPZGeoMesh *gmesh, int nDiv);
template <class TGeo>
void CompareSameDimensionNonLinNeighbour(int nDiv);
}

TEST_CASE("geoblend_tests","[blend_tests]") {
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    {
        const int nDiv = 4;
        auto gmesh = TPZGenSpecialGrid::CreateGeoMesh2D_Circle(0);
        blendtest::TestMesh2D(gmesh,nDiv);
        delete gmesh;
    }
    
    gRefDBase.InitializeUniformRefPattern(ETetraedro);
    gRefDBase.InitializeUniformRefPattern(EPiramide);
    gRefDBase.InitializeUniformRefPattern(EPrisma);
    gRefDBase.InitializeUniformRefPattern(ECube);
    {
        const int nDiv = 2;
        auto gmesh = TPZGenSpecialGrid::CreateGeoMesh3D_DividedSphere(0);
        blendtest::TestMesh3D(gmesh, nDiv);
        delete gmesh;
    }
    {
        const int nDiv = 2;
        blendtest::CompareSameDimensionNonLinNeighbour<pzgeom::TPZGeoLinear>(nDiv);
        blendtest::CompareSameDimensionNonLinNeighbour<pzgeom::TPZGeoTriangle>(nDiv);
        blendtest::CompareSameDimensionNonLinNeighbour<pzgeom::TPZGeoQuad>(nDiv);
    }
}


TEST_CASE("compare_blend_quad","[blend_tests]") {
    blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoQuad>();
    blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoTriangle>();
    blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoTetrahedra>();
    blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoCube>();
    blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoPrism>();
    blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoPyramid>();
}

TEST_CASE("quadrilateral_blend_semicircle","[blend_tests]") {
    //defining the analytical solution against which the mapping will be compared. the quadrilateral will be mapped
    //to a quarter of a ring. the inner and outer radii are, respectively, 0.5 and 1.0.
    const REAL innerRadius = 0.5;
    const REAL outerRadius = 1.0;
    auto analyticX = [innerRadius, outerRadius ] (const TPZVec<REAL> &pt) {
        TPZVec<REAL> map(3,-1);
        map[0] = innerRadius * (1. - pt[0]) * cos((1. + pt[1]) * M_PI_4) * 0.5 +
        outerRadius * (1. + pt[0]) * cos((1. + pt[1]) * M_PI_4) * 0.5;
        map[1] = innerRadius * (1. - pt[0]) * sin((1. + pt[1]) * M_PI_4) * 0.5 +
        outerRadius * (1. + pt[0]) * sin((1. + pt[1]) * M_PI_4) * 0.5;
        map[2] = 0.0;
        return map;
    };
    
    auto analyticGradX = [innerRadius, outerRadius ] (const TPZVec<REAL> &pt) {
        TPZFMatrix<REAL> grad(3,2,-1);
        grad(0,0) = -0.5 * innerRadius * cos((1. + pt[1]) * M_PI_4)
        + 0.5 * outerRadius * cos((1. + pt[1]) * M_PI_4);
        grad(0,1) = -0.5 * M_PI_4 * innerRadius * (1. - pt[0]) * sin((1. + pt[1]) * M_PI_4)
        - 0.5 * M_PI_4 * outerRadius * (1. + pt[0]) * sin((1. + pt[1]) * M_PI_4);
        grad(1,0) = -0.5 * innerRadius * sin((1. + pt[1]) * M_PI_4)
        + 0.5 * outerRadius * sin((1. + pt[1]) * M_PI_4);
        grad(1,1) =  0.5 * M_PI_4 * innerRadius * (1. - pt[0]) * cos((1. + pt[1]) * M_PI_4)
        + 0.5 * M_PI_4 * outerRadius * (1. + pt[0]) * cos((1. + pt[1]) * M_PI_4);
        grad(2,0) = 0.0;
        grad(2,1) = 0.0;
        return grad;
    };
    
    TPZGeoMesh *gmesh = TPZGenSpecialGrid::CreateArcMesh(innerRadius,outerRadius);

    auto quad = gmesh->Element(0);
    const int64_t nPoints = 20;
    const REAL tol = blendtest::tol;
    
    TPZManVector<REAL,2> qsi(2);
    TPZManVector<REAL,3> xBlend(3);
    TPZManVector<REAL,3> xAnalytic(3);
    TPZFNMatrix<9,REAL> gradxBlend(3,3);
    TPZFNMatrix<9,REAL> gradxAnalytic(3,3);
    bool hasAnErrorOccurred = false;
    int nPointsTotal = 0, nErrorsX = 0, nErrorsGradX = 0;
    for (int iPt = 0; iPt < nPoints; ++iPt) {
        qsi[0] = -1. + 2 * (((REAL) iPt)/(nPoints-1));
        for(int jPt = 0; jPt < nPoints; ++jPt) {
            nPointsTotal++;
            qsi[1] = -1. + 2 * (((REAL) jPt) / (nPoints - 1));
            quad->X(qsi, xBlend);
            xAnalytic = analyticX(qsi);
            hasAnErrorOccurred =
            blendtest::CheckVectors(xBlend, "xBlend", xAnalytic, "xAnalytic", blendtest::tol);
            REQUIRE(!hasAnErrorOccurred);
            if(hasAnErrorOccurred) nErrorsX++;
            quad->GradX(qsi, gradxBlend);
            gradxAnalytic = analyticGradX(qsi);
            hasAnErrorOccurred =
            blendtest::CheckMatrices(gradxBlend, "gradxBlend", gradxAnalytic, "gradxAnalytic", blendtest::tol);
            REQUIRE(!hasAnErrorOccurred);
            if(hasAnErrorOccurred) nErrorsGradX ++;
            //TODO: (IMPORTANT) deactivate the log after testing
        }
    }
#ifdef BLEND_VERBOSE
    std::cout<<"quadrilateral_blend_semicircle:"<<'\n';
    std::cout<<"\tnPointsTotal = "<<nPointsTotal<<'\n';
    std::cout<<"\tnErrorsX = "<<nErrorsX<<'\n';
    std::cout<<"\tnErrorsGradX = "<<nErrorsGradX<<std::endl;
#endif
    delete gmesh;
}
namespace blendtest{

bool CheckVectors(const TPZVec<REAL> &x1, std::string name1,
                  const TPZVec<REAL> &x2, std::string name2, const REAL &tol){
    const auto VAL_WIDTH = 15;
    REAL diff = 0;
    bool hasAnErrorOccurred = false;
    if(x1.size() != x2.size()){
        hasAnErrorOccurred = true;
        return hasAnErrorOccurred;
    }
    for (int i = 0; i < x1.size(); i++) {
        diff += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    if (diff > tol) {
        hasAnErrorOccurred = true;
    }
    if (hasAnErrorOccurred) {
#ifdef BLEND_VERBOSE
        std::ostringstream x1m, x2m;
        x1m << name1 << '\n';
        x2m << name2 << '\n';
        for (int i = 0; i < x1.size(); i++) {
            x1m << std::setw(VAL_WIDTH) << std::right << x1[i] << "\t";
            x2m << std::setw(VAL_WIDTH) << std::right << x2[i] << "\t";
        }
        std::cout << x1m.str() << '\n';
        std::cout << x2m.str() << '\n';
        std::cout << "diff :" << diff << std::endl;
#endif
    }
    return hasAnErrorOccurred;
}

bool CheckMatrices(const TPZFMatrix<REAL> &gradx1, std::string name1,
                   const TPZFMatrix<REAL> &gradx2, std::string name2, const REAL &tol){
    const auto VAL_WIDTH = 15;
    REAL diff = 0;
    bool hasAnErrorOccurred = false;
    if((gradx1.Rows() != gradx2.Rows()) || (gradx1.Cols() != gradx2.Cols())){
        hasAnErrorOccurred = true;
        return hasAnErrorOccurred;
    }
    for (int i = 0; i < gradx1.Rows(); i++) {
        for (int j = 0; j < gradx1.Cols(); j++) {
            diff += (gradx1.GetVal(i, j) - gradx2.GetVal(i, j)) * (gradx1.GetVal(i, j) - gradx2.GetVal(i, j));
        }
    }
    if (diff > tol) {
        hasAnErrorOccurred = true;
    }
    if (hasAnErrorOccurred) {
#ifdef BLEND_VERBOSE
        std::ostringstream x1m, x2m;
        x1m << name1 << std::endl;
        x2m << name2 << std::endl;
        for (int i = 0; i < gradx1.Rows(); i++) {
            for (int j = 0; j < gradx1.Cols(); j++) {
                x1m << std::setw(VAL_WIDTH) << std::right
                << gradx1.GetVal(i, j) << "\t";
                x2m << std::setw(VAL_WIDTH) << std::right
                << gradx2.GetVal(i, j) << "\t";
            }
            x1m << '\n';
            x2m << '\n';
        }
        std::cout << x1m.str() << '\n';
        std::cout << x2m.str() << '\n';
        std::cout << "diff :" << diff << std::endl;
#endif
    }
    return hasAnErrorOccurred;
}

template <class TGeo>
void CompareQuadraticAndBlendEls() {
#ifdef BLEND_VERBOSE
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
    auto elType = TGeo::Type();
    
    const int64_t nCornerNodes = TGeo::NCornerNodes;
    const int64_t nEdges = TGeo::NumSides(1);
    const int64_t dim = TGeo::Dimension;
    const REAL sphereRadius = 1;
    
    TPZVec<REAL> phiPts(nCornerNodes,-1),thetaPts(nCornerNodes,-1); //r is always equal to sphereRadius
    switch(elType){
        case ETriangle:
            for(int i = 0; i < nCornerNodes; i++){
                thetaPts[i] = M_PI/2;
                phiPts[i] = i * M_PI/3;
            }
            break;
        case EQuadrilateral:
            for(int i = 0; i < nCornerNodes; i++){
                thetaPts[i] = M_PI/2;
                phiPts[i] = i * M_PI/2;
            }
            break;
        case ETetraedro:
            for(int i = 0; i < 3; i++){
                thetaPts[i] = 3 * M_PI/4;
                phiPts[i] = i * M_PI/3;
            }
            thetaPts[3] = 0;
            phiPts[3] = 0;
            break;
        case EPrisma:
            for(int i = 0; i < 3; i++){
                thetaPts[i] = 3 * M_PI/4;
                phiPts[i] = i * M_PI/3;
            }
            for(int i = 0; i < 3; i++){
                thetaPts[i+3] = 1 * M_PI/4;
                phiPts[i+3] = i * M_PI/3;
            }
            break;
        case ECube:
            for(int i = 0; i < 4; i++){
                thetaPts[i] = 3 * M_PI/4;
                phiPts[i] = i * M_PI/2;
            }
            for(int i = 0; i < 4; i++){
                thetaPts[4+i] = 1 * M_PI/4;
                phiPts[4+i] = i * M_PI/2;
            }
            break;
        case EPiramide:
            for(int i = 0; i < 4; i++){
                thetaPts[i] = 3 * M_PI/4;
                phiPts[i] = i * M_PI/2;
            }
            thetaPts[4] = 0;
            phiPts[4] = 0;
            break;
        default:
            DebugStop();
    }
    
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    TPZManVector<REAL,3> coordsOffset(3,0);
    coordsOffset[0] = 2;
    TPZManVector<int, 12> midSideNodesIndexVec(nEdges, -1);
    const int nNodesOriginalMesh = nCornerNodes + nEdges;
    {
        TPZManVector<REAL, 3> coord(3, 0.);
        
        ///CREATE NODES FOR QUADRATIC ELEMENT
        for (int64_t i = 0; i < nCornerNodes; i++) {
            const REAL theta = ((REAL) i / nCornerNodes) * 2 * M_PI;
            coord[0] = sphereRadius * sin(thetaPts[i]) * cos(phiPts[i]);
            coord[1] = sphereRadius * sin(thetaPts[i]) * sin(phiPts[i]);
            coord[2] = sphereRadius * cos(thetaPts[i]);
            //                std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
        }//quadraticElCornerNodes
        
        
        //CREATE MIDSIDE NODES
        TPZStack<int> edgeStack;
        TGeo::LowerDimensionSides(TGeo::NSides - 1, edgeStack, 1);
        for (int64_t edgeIndex = 0;
             edgeIndex < edgeStack.NElements(); edgeIndex++) {
            const int64_t edge = edgeStack[edgeIndex];
            const int64_t nNodesSide = TGeo::NSideNodes(edge);
            coord[0] = 0;
            coord[1] = 0;
            coord[2] = 0;
            for (int64_t node = 0; node < nNodesSide; node++) {
                const int64_t nodeIndex = TGeo::SideNodeLocId(edge, node);
                coord[0] += gmesh->NodeVec()[nodeIndex].Coord(0);
                coord[1] += gmesh->NodeVec()[nodeIndex].Coord(1);
                coord[2] += gmesh->NodeVec()[nodeIndex].Coord(2);
            }
            const REAL norm = sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
            coord[0] *= sphereRadius / norm;
            coord[1] *= sphereRadius / norm;
            coord[2] *= sphereRadius / norm;
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
            midSideNodesIndexVec[edgeIndex] = newindex;
        }
        
        ///CREATE NODES FOR BLEND ELEMENT
        for (int i = 0; i < nNodesOriginalMesh; i++) {
            coord[0] = gmesh->NodeVec()[i].Coord(0) + coordsOffset[0];
            coord[1] = gmesh->NodeVec()[i].Coord(1) + coordsOffset[1];
            coord[2] = gmesh->NodeVec()[i].Coord(2) + coordsOffset[2];
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
        }
    }
    
    const int matIdVol = 1, matIdSphere = 2;
    TPZManVector<int64_t,8> nodesIdVec(1);
    int64_t elId = 0;
    
    ///////CREATE MAX DIM ELEMENTS
    nodesIdVec.Resize(nCornerNodes + nEdges);
    for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i;
    
    TPZGeoEl *quadraticEl = nullptr;
    switch(elType){
        case ETriangle:
            quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(elId,nodesIdVec,matIdVol,*gmesh);
            break;
        case EQuadrilateral:
            quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(elId,nodesIdVec,matIdVol,*gmesh);
            break;
        case ETetraedro:
            quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticTetra>(elId,nodesIdVec,matIdVol,*gmesh);
            break;
        case EPrisma:
            quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticPrism>(elId,nodesIdVec,matIdVol,*gmesh);
            break;
        case ECube:
            quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticCube>(elId,nodesIdVec,matIdVol,*gmesh);
            break;
        case EPiramide:
            quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticPyramid>(elId,nodesIdVec,matIdVol,*gmesh);
            break;
        default:
            DebugStop();
            break;
    }
    nodesIdVec.Resize(nCornerNodes);
    for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i + nNodesOriginalMesh;
    TPZGeoEl *blendEl = new TPZGeoElRefLess<pzgeom::TPZGeoBlend<TGeo>>(elId,nodesIdVec,matIdVol,*gmesh);
    
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "gmesh_" + TGeo::TypeName() + ".txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
    ///////CREATE MAX DIM - 1 ELEMENTS FOR BLENDING
    {
        TPZStack<int> containedSides;
        TGeo::LowerDimensionSides(TGeo::NSides - 1, containedSides, TGeo::Dimension - 1);
        for (int64_t subSideIndex = 0;
             subSideIndex < containedSides.NElements(); subSideIndex++){
            const int64_t subSide = containedSides[subSideIndex];
            const int64_t nNodesSide = TGeo::NSideNodes(subSide);
            TPZStack<int> containedEdges;
            TGeo::LowerDimensionSides(subSide, containedEdges, 1);
            if(TGeo::SideDimension(subSide) == 1)   containedEdges.Push(subSide);
            const int64_t nEdgesSide = containedEdges.NElements();
            nodesIdVec.Resize(nNodesSide+nEdgesSide);
            TPZManVector<int64_t,8> nodesIdVecQuad(nNodesSide+nEdgesSide);
            for (int64_t node = 0; node < nNodesSide; node++){
                const int64_t nodeId = TGeo::SideNodeLocId(subSide,node);
                nodesIdVec[node] = blendEl->Node(nodeId).Id();
                nodesIdVecQuad[node] = quadraticEl->Node(nodeId).Id();
                //                    std::cout<<nodesIdVec[node]<<std::endl;
            }
            for (int64_t edge = 0; edge < nEdgesSide; edge++){
                const int64_t edgeId = containedEdges[edge];
                nodesIdVec[nNodesSide + edge] = nNodesOriginalMesh + midSideNodesIndexVec[edgeId- nCornerNodes];
                nodesIdVecQuad[nNodesSide + edge] = midSideNodesIndexVec[edgeId- nCornerNodes];
                //                    std::cout<<nodesIdVec[nNodesSide + edge]<<std::endl;
            }
            switch (TGeo::Type(subSide)){
                case EOned:
                    new TPZGeoElRefLess<pzgeom::TPZQuadraticLine>(nodesIdVec,matIdSphere,*gmesh);
                    new TPZGeoElRefLess<pzgeom::TPZQuadraticLine>(nodesIdVecQuad,matIdSphere,*gmesh);
                    break;
                case ETriangle:
                    new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(nodesIdVec,matIdSphere,*gmesh);
                    new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(nodesIdVecQuad,matIdSphere,*gmesh);
                    break;
                case EQuadrilateral:
                    new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(nodesIdVec,matIdSphere,*gmesh);
                    new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(nodesIdVecQuad,matIdSphere,*gmesh);
                    break;
                default:
                    DebugStop();
            }
        }
    }
    gmesh->BuildConnectivity();
    
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "gmesh_" + TGeo::TypeName() + ".txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "gmesh_" + TGeo::TypeName() + ".vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
    //X COMPARE
    {
        TPZManVector<REAL,3> xi;
        REAL notUsedHere = -1;
        uint64_t errorsInterior = 0;
        TPZManVector<uint64_t,18> errorsSide(TGeo::NSides - TGeo::NNodes - 1,0);
        uint64_t nPoints = 0;
        for(int iSide = TGeo::NNodes; iSide < TGeo::NSides; iSide ++){
            auto intRule = blendEl->CreateSideIntegrationRule(iSide, pOrder);
            xi.Resize(dim,0);
            for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                bool hasAnErrorOccurred = false;
                nPoints++;
                TPZManVector<REAL,3> xiSide(TGeo::SideDimension(iSide),0);
                intRule->Point(iPt,xiSide,notUsedHere);
                auto transf = TGeo::TransformSideToElement(iSide);
                transf.Apply(xiSide,xi);
                TPZManVector<REAL,3> xBlend(3);
                blendEl->X(xi, xBlend);
                for(int iX = 0; iX < xBlend.size(); iX++) xBlend[iX] -= coordsOffset[iX];
                TPZManVector<REAL,3> xQuad(3);
                quadraticEl->X(xi, xQuad);
                hasAnErrorOccurred = blendtest::CheckVectors(xBlend,"xBlend",xQuad,"xQuad",blendtest::tol);
                REQUIRE(!hasAnErrorOccurred);
                if(hasAnErrorOccurred){
#ifdef BLEND_VERBOSE
                    if(iSide < TGeo::NSides - 1) errorsSide[iSide - TGeo::NNodes]+=1;
                    if(iSide == TGeo::NSides - 1) errorsInterior++;
                    const auto VAL_WIDTH = 15;
                    if(iSide < TGeo::NSides - 1){
                        auto neigh = blendEl->Neighbour(iSide);
                        if(neigh.Id() != blendEl->Id()) {
                            TPZGeoElSide thisside(blendEl, iSide);
                            auto neighTransf = thisside.NeighbourSideTransform(neigh);
                            TPZManVector<REAL, 3> neighXi(neigh.Dimension(), 0);
                            neighTransf.Apply(xiSide, neighXi);
                            TPZManVector<REAL, 3> xNeigh(3);
                            neigh.X(neighXi, xNeigh);
                            std::cout<<"x_neigh_blend:"<<std::endl;
                            
                            for (int i = 0; i < xNeigh.size(); i++) {
                                std::cout << std::setw(VAL_WIDTH) << std::right << xNeigh[i] - coordsOffset[i]
                                << "\t";
                            }
                            std::cout<<std::endl;
                        }
                    }
                    
                    if(iSide < TGeo::NSides - 1){
                        auto neigh = quadraticEl->Neighbour(iSide);
                        if(neigh.Id() != quadraticEl->Id()) {
                            
                            TPZGeoElSide thisside(quadraticEl, iSide);
                            auto neighTransf = thisside.NeighbourSideTransform(neigh);
                            TPZManVector<REAL, 3> neighXi(neigh.Dimension(), 0);
                            neighTransf.Apply(xiSide, neighXi);
                            TPZManVector<REAL, 3> xNeigh(3);
                            neigh.X(neighXi, xNeigh);
                            
                            std::cout<<"x_neigh_quad:"<<std::endl;
                            for (int i = 0; i < xNeigh.size(); i++) {
                                std::cout << std::setw(VAL_WIDTH) << std::right << xNeigh[i]<< "\t";
                            }
                            std::cout<<std::endl;
                        }
                    }
#endif
                }
            }
#ifdef BLEND_VERBOSE
            if(iSide == TGeo::NSides - 1){
                std::cout<<"\tNumber of interior points: "<<intRule->NPoints()<<"\tErrors: "<<errorsInterior<<std::endl;
            }else{
                std::cout<<"\tNumber of side points for side "<<iSide<<": "<<intRule->NPoints()<<"\tErrors: "<<errorsSide[iSide - TGeo::NNodes]<<std::endl;
            }
#endif
        }
        uint64_t errorsTotal = errorsInterior;
        for(int i = 0; i< errorsSide.size(); i++){
            errorsTotal +=errorsSide[i];
        }
#ifdef BLEND_VERBOSE
        std::cout<<"\tNumber of points: "<<nPoints<<"\tErrors: "<<errorsTotal<<std::endl;
#endif
    }
    delete gmesh;
}


template <class TGeo>
void CompareSameDimensionNonLinNeighbour(int nref) {
#ifdef BLEND_VERBOSE
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
    auto elType = TGeo::Type();
    
    const int64_t nCornerNodes = TGeo::NCornerNodes;
    const int64_t nEdges = TGeo::NumSides(1);
    const int64_t dim = TGeo::Dimension;
    const REAL sphereRadius = 1;
    TPZManVector<REAL,3> sphereCenter(3,0);
    TPZVec<REAL> phiPts(nCornerNodes,-1),thetaPts(nCornerNodes,-1); //r is always equal to sphereRadius
    
    switch(elType){
        case EOned:
            for(int i = 0; i < nCornerNodes; i++){
                thetaPts[i] = M_PI/2;
                phiPts[i] = i * M_PI/2;
            }
            break;
        case ETriangle:
            for(int i = 0; i < nCornerNodes; i++){
                thetaPts[i] = M_PI/2;
                phiPts[i] = i * M_PI/2;
            }
            break;
        case EQuadrilateral:
            for(int i = 0; i < nCornerNodes; i++){
                thetaPts[i] = M_PI/2;
                phiPts[i] = i * M_PI/2;
            }
            break;
        default:
            DebugStop();
    }
    
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    TPZManVector<int, 12> midSideNodesIndexVec(nEdges, -1);
    const int nNodesOriginalMesh = nCornerNodes + nEdges;
    {
        TPZManVector<REAL, 3> coord(3, 0.);
        
        ///CREATE NODES FOR LINEAR ELEMENT
        for (int64_t i = 0; i < nCornerNodes; i++) {
            const REAL theta = ((REAL) i / nCornerNodes) * 2 * M_PI;
            coord[0] = sphereRadius * sin(thetaPts[i]) * cos(phiPts[i]);
            coord[1] = sphereRadius * sin(thetaPts[i]) * sin(phiPts[i]);
            coord[2] = sphereRadius * cos(thetaPts[i]);
            //                std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
        }//
        
        
        //CREATE MIDSIDE NODES
        TPZStack<int,TGeo::NSides> edgeStack;
        for(int iSide = 0; iSide < TGeo::NSides; iSide++){
            if(TGeo::SideDimension(iSide) == 1) edgeStack.Push(iSide);
        }
        //            TGeo::LowerDimensionSides(TGeo::NSides - 1, edgeStack, 1);
        for (int64_t edgeIndex = 0;
             edgeIndex < edgeStack.NElements(); edgeIndex++) {
            const int64_t edge = edgeStack[edgeIndex];
            const int64_t nNodesSide = TGeo::NSideNodes(edge);
            REAL sumPhiNodes = 0;
            REAL sumThetaNodes = 0;
            for (int64_t node = 0; node < nNodesSide; node++) {
                const int64_t nodeIndex = TGeo::SideNodeLocId(edge, node);
                sumPhiNodes += phiPts[nodeIndex];
                sumThetaNodes += thetaPts[nodeIndex];
            }
            sumPhiNodes /= nNodesSide;
            sumThetaNodes /= nNodesSide;
            coord[0] = sphereRadius * sin(sumThetaNodes) * cos(sumPhiNodes);
            coord[1] = sphereRadius * sin(sumThetaNodes) * sin(sumPhiNodes);
            coord[2] = sphereRadius * cos(sumThetaNodes);
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
            midSideNodesIndexVec[edgeIndex] = newindex;
        }
    }
    
    const int matIdVol = 1, matIdSphere = 2;
    TPZManVector<int64_t,8> nodesIdVec(1);
    int64_t elId = 0;
    
    ///////CREATE MAX DIM ELEMENTS
    nodesIdVec.Resize(nNodesOriginalMesh);
    for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i;
    
    TPZGeoEl *nonLinearEl = nullptr;
    switch(elType){
        case EOned:
            nonLinearEl = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(elId,nodesIdVec,matIdVol, *gmesh);
            break;
        case ETriangle:
            nonLinearEl = new TPZGeoElRefPattern<pzgeom::TPZTriangleSphere<pzgeom::TPZGeoTriangle>>(elId,nodesIdVec,matIdVol,*gmesh);
        {
            auto sphere = dynamic_cast<TPZGeoElRefPattern<pzgeom::TPZTriangleSphere<pzgeom::TPZGeoTriangle>> *> (nonLinearEl);
            sphere->Geom().SetData(sphereRadius, sphereCenter);
        }
            break;
        case EQuadrilateral:
            nonLinearEl = new TPZGeoElRefPattern<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad>>(elId,nodesIdVec,matIdVol,*gmesh);
        {
            auto sphere = dynamic_cast<TPZGeoElRefPattern<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad>> *> (nonLinearEl);
            sphere->Geom().SetData(sphereRadius, sphereCenter);
        }
            break;
        default:
            DebugStop();
            break;
    }
    nodesIdVec.Resize(nCornerNodes);
    for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i;
    TPZGeoEl *blendEl = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<TGeo>>(elId,nodesIdVec,matIdVol,*gmesh);
    
    gmesh->BuildConnectivity();
    {
        TPZVec<TPZGeoEl *> sons;
        for(int iDiv = 0; iDiv < nref; iDiv++){
            const int nel = gmesh-> NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh-> Element(iel);
                if (geo->NSubElements() == 0) {
                    geo->Divide(sons);
                }
            }
        }
    }
    
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "gmesh_nonlin_" + TGeo::TypeName() + ".txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "gmesh_nonlin_" + TGeo::TypeName()+ ".vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
    
    //X COMPARE
    {
        TPZManVector<REAL,3> xi;
        REAL notUsedHere = -1;
        uint64_t errorsInterior = 0;
        TPZManVector<uint64_t,18> errorsSide(TGeo::NSides - TGeo::NNodes - 1,0);
        uint64_t nPoints = 0;
        for(int iSide = TGeo::NNodes; iSide < TGeo::NSides; iSide ++){
            auto intRule = blendEl->CreateSideIntegrationRule(iSide, pOrder);
            xi.Resize(dim,0);
            for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                bool hasAnErrorOccurred = false;
                nPoints++;
                TPZManVector<REAL,3> xiSide(TGeo::SideDimension(iSide),0);
                intRule->Point(iPt,xiSide,notUsedHere);
                auto transf = TGeo::TransformSideToElement(iSide);
                transf.Apply(xiSide,xi);
                TPZManVector<REAL,3> xBlend(3);
                blendEl->X(xi, xBlend);
                TPZManVector<REAL,3> xNonLin(3);
                nonLinearEl->X(xi, xNonLin);
                hasAnErrorOccurred = blendtest::CheckVectors(xBlend, "xBlend", xNonLin, "xNonLin", blendtest::tol);
                REQUIRE(!hasAnErrorOccurred);
                if(hasAnErrorOccurred){
#ifdef BLEND_VERBOSE
                    if(iSide < TGeo::NSides - 1) errorsSide[iSide - TGeo::NNodes]+=1;
                    if(iSide == TGeo::NSides - 1) errorsInterior++;
                    const auto VAL_WIDTH = 15;
                    if(iSide < TGeo::NSides - 1){
                        auto neigh = blendEl->Neighbour(iSide);
                        if(neigh.Id() != blendEl->Id()) {
                            TPZGeoElSide thisside(blendEl, iSide);
                            auto neighTransf = thisside.NeighbourSideTransform(neigh);
                            TPZManVector<REAL, 3> neighXi(neigh.Dimension(), 0);
                            neighTransf.Apply(xiSide, neighXi);
                            TPZManVector<REAL, 3> xNeigh(3);
                            neigh.X(neighXi, xNeigh);
                            std::cout<<"x_neigh_blend:"<<std::endl;
                            for (int i = 0; i < xNeigh.size(); i++) {
                                std::cout << std::setw(VAL_WIDTH) << std::right << xNeigh[i]
                                << "\t";
                            }
                            std::cout<<std::endl;
                        }
                    }
                    
                    if(iSide < TGeo::NSides - 1){
                        auto neigh = nonLinearEl->Neighbour(iSide);
                        if(neigh.Id() != nonLinearEl->Id()) {
                            TPZGeoElSide thisside(nonLinearEl, iSide);
                            auto neighTransf = thisside.NeighbourSideTransform(neigh);
                            TPZManVector<REAL, 3> neighXi(neigh.Dimension(), 0);
                            neighTransf.Apply(xiSide, neighXi);
                            TPZManVector<REAL, 3> xNeigh(3);
                            neigh.X(neighXi, xNeigh);
                            std::cout<<"x_neigh_quad:"<<std::endl;
                            for (int i = 0; i < xNeigh.size(); i++) {
                                std::cout << std::setw(VAL_WIDTH) << std::right << xNeigh[i]<< "\t";
                            }
                            std::cout<<std::endl;
                        }
                    }
#endif
                }
            }
#ifdef BLEND_VERBOSE
            if(iSide == TGeo::NSides - 1){
                std::cout<<"\tNumber of interior points: "<<intRule->NPoints()<<"\tErrors: "<<errorsInterior<<std::endl;
            }else{
                std::cout<<"\tNumber of side points for side "<<iSide<<": "<<intRule->NPoints()<<"\tErrors: "<<errorsSide[iSide - TGeo::NNodes]<<std::endl;
            }
#endif
        }
        uint64_t errorsTotal = errorsInterior;
        for(int i = 0; i< errorsSide.size(); i++){
            errorsTotal +=errorsSide[i];
        }
#ifdef BLEND_VERBOSE
        std::cout<<"\tNumber of points: "<<nPoints<<"\tErrors: "<<errorsTotal<<std::endl;
#endif
    }
    delete gmesh;
}

void TestMesh2D(TPZGeoMesh *gmesh, int nDiv)
{
    {
        TPZManVector<REAL,3> xiReal;
        TPZManVector<Fad<REAL>,3> xiFad;
        REAL weight = -1;//useless
        const int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->ElementVec()[iel];
            if (geo && !geo->HasSubElement()) {
                uint64_t errors = 0;
                auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                xiReal.Resize(geo->Dimension(),0);
                xiFad.Resize(geo->Dimension(),0);
                for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                    bool hasAnErrorOccurred = false;
                    intRule->Point(iPt,xiReal,weight);
                    for(int x = 0; x < geo->Dimension(); x++){
                        xiFad[x] = Fad<REAL>(geo->Dimension(),x,xiReal[x]);
                    }
                    //                        Fad<REAL> func;
                    //                        func = 9*xiFad[0];
                    //                        for(int i =0; i < func.size(); i++){
                    //                            std::cout<<"dx["<<i<<"]:\t"<<func.dx(i)<<std::endl;
                    //                        }
                    TPZManVector<Fad<REAL>,3> xFad(3);
                    geo->X(xiFad, xFad);
                    TPZFNMatrix<9,REAL> gradXreal(3,3,0);
                    geo->GradX(xiReal, gradXreal);
                    TPZFNMatrix<9,REAL> gradXfad(gradXreal.Rows(),gradXreal.Cols(),0.0);
                    for(int i = 0; i < gradXreal.Rows(); i++){
                        for(int j = 0; j < gradXreal.Cols(); j++){
                            gradXfad(i,j) = xFad[i].dx(j);
                        }
                    }
                    hasAnErrorOccurred = blendtest::CheckMatrices(gradXreal,"gradXreal",gradXfad,"gradXfad",blendtest::tol);
                    REQUIRE(!hasAnErrorOccurred);
                    if(hasAnErrorOccurred){
                        errors++;
                    }
                }
#ifdef BLEND_VERBOSE
                if(errors > 0 || geo->IsGeoBlendEl()){
                    std::cout << "============================"
                    << std::endl;
                    std::cout
                    << "Element: " << geo->Id()
                    << "\tType: " << MElementType_Name(geo->Type());
                    std::cout << "\tIs blend? : " << geo->IsGeoBlendEl()
                    << std::endl;
                    std::cout
                    << "\tNumber of points: " << intRule->NPoints()
                    << "\tErrors: " << errors << std::endl;
                }
#endif
            }
        }
    }
    
    
    {
        TPZVec<TPZGeoEl *> sons;
        std::vector<std::string> loading = {"-","/","|","\\"};
        for (int iDiv = 0; iDiv < nDiv; iDiv++) {
            std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                std::cout<<"\b"<<loading[iel%4]<<std::flush;
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    geo->Divide(sons);
                }
            }
            std::cout<<"\b";
        }
    }
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh2D.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh2D.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
}


void TestMesh3D(TPZGeoMesh *gmesh, int nDiv)
{
    {
        TPZManVector<REAL,3> xiReal;
        TPZManVector<Fad<REAL>,3> xiFad;
        REAL weight = -1;//useless
        const int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->ElementVec()[iel];
            if (geo && !geo->HasSubElement()) {
                uint64_t errors = 0;
                auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                xiReal.Resize(geo->Dimension(),0);
                xiFad.Resize(geo->Dimension(),0);
                for(int iPt = 0; iPt < intRule->NPoints(); iPt++) {
                    bool hasAnErrorOccurred = false;
                    intRule->Point(iPt, xiReal, weight);
                    for (int x = 0; x < geo->Dimension(); x++) {
                        xiFad[x] = Fad<REAL>(geo->Dimension(), x, xiReal[x]);
                    }
                    
                    TPZManVector<Fad<REAL>, 3> xFad(3);
                    geo->X(xiFad, xFad);
                    TPZFNMatrix<9, REAL> gradXreal(3, 3, 0);
                    geo->GradX(xiReal, gradXreal);
                    TPZFNMatrix<9,REAL> gradXfad(gradXreal.Rows(),gradXreal.Cols(),0.0);
                    for (int i = 0; i < gradXreal.Rows(); i++) {
                        for (int j = 0; j < gradXreal.Cols(); j++) {
                            gradXfad(i, j) = xFad[i].dx(j);
                        }
                    }
                    hasAnErrorOccurred = blendtest::CheckMatrices(gradXreal, "gradXreal", gradXfad, "gradXfad",
                                                                  blendtest::tol);
                    REQUIRE(!hasAnErrorOccurred);
                    if (hasAnErrorOccurred) {
                        errors++;
                    }
                }
#ifdef BLEND_VERBOSE
                if(errors > 0 || geo->IsGeoBlendEl()){
                    std::cout << "============================"
                    << std::endl;
                    std::cout
                    << "Element: " << geo->Id()
                    << "\tType: " << MElementType_Name(geo->Type());
                    std::cout << "\tIs blend? : " << geo->IsGeoBlendEl()
                    << std::endl;
                    std::cout
                    << "\tNumber of points: " << intRule->NPoints()
                    << "\tErrors: " << errors << std::endl;
                }
#endif
            }
        }
    }
    {
        TPZVec<TPZGeoEl *> sons;
        std::vector<std::string> loading = {"-","/","|","\\"};
        for (int iDiv = 0; iDiv < nDiv; iDiv++) {
#ifdef BLEND_OUTPUT_TXT
            {
                const std::string meshFileName =
                "blendmesh3D_ref" + std::to_string(iDiv) + ".txt";
                std::ofstream outTXT(meshFileName.c_str());
                gmesh->Print(outTXT);
                outTXT.close();
            }
#endif
            std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                std::cout<<"\b"<<loading[iel%4]<<std::flush;
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    geo->Divide(sons);
                }
            }
            std::cout<<"\b";
        }
    }
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh3D.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh3D.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
}


}
