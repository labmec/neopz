/**
 * @file HCurlUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of HCurl-conforming approximation spaces
 *
 */
#include <iostream>

#include "TPZTopologyUtils.h"
#include "TPZGeoMeshTools.h"

#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

#include "pzgeoelrefless.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testhcurl");
#endif

#include "fad.h"

#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "pzcmesh.h"
#include "Projection/TPZHCurlProjection.cpp"
#include "TPZLinearAnalysis.h"
#include "pzintel.h"
#include "TPZCompElHCurl.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZVTKGeoMesh.h"


#include <catch2/catch.hpp>

#define REQUIRE_MESSAGE(cond, msg) do { INFO(msg); REQUIRE(cond); } while((void)0, 0)

//#define VERBOSE_HCURL //outputs useful debug info
//#define HCURL_OUTPUT_TXT//export mesh in .txt file (for debugging)
//#define HCURL_OUTPUT_VTK//export mesh in .vtk file (for debugging)

namespace hcurltest{
    constexpr REAL tol = 1e-10;
    /**
     * This unit test aims to verify the vectors used to build the HCurl approximation space.
     * It was rewritten based on TestFunctionTracesUniformMesh, therefore it could be simplified.
     * @param type element type (triangle, quadrilateral, etc.
     * @param dim dimension of the element
     */
    void TestVectorTracesUniformMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                                     MMeshType type);
    /**
     * This unit test aims to verify the trace compatibility of the HCurl approximation space in a UNIFORM mesh.
     * @param type element type (triangle, quadrilateral, etc.
     * @param dim dimension of the element
     */
    void TestFunctionTracesUniformMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                                       MMeshType type, const int pOrder);
    void TestFunctionCurlUniformMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                                       MMeshType type, const int pOrder);
    //auxiliary funcs
    template <class TGEOM>
    void PrintShapeFunctions(const int pOrder);

    TPZAutoPointer<TPZCompMesh> CreateCMesh(MMeshType type, const int pOrder);
    
    bool CheckTracesFunc(REAL &diffTrace, const TPZVec<REAL> &elShapeFunc, const TPZVec<REAL>&neighShapeFunc,
                         const TPZVec<REAL> &vec, const int sideDim, TPZVec<REAL> &elTrace, TPZVec<REAL> &neighTrace);
    
    void VectorProduct(const TPZVec<REAL> &, const TPZVec<REAL> &, TPZVec<REAL> &);
}


TEST_CASE("Testing trace of HCurl functions",
          "[hcurl_tests][mesh][topology]") {
    
    // SECTION("PrintFunctions"){
    //     constexpr int pOrder{3};
    //     hcurltest::PrintShapeFunctions<pztopology::TPZQuadrilateral>(pOrder);
    // }
    
    constexpr int pOrder{1};
    constexpr int maxK{5};
    auto meshType = GENERATE(MMeshType::ETriangular,
                             MMeshType::EQuadrilateral,
                             MMeshType::ETetrahedral,
                             MMeshType::EHexahedral
                             // MMeshType::EPrismatic//NEEDSFIX
                             );
    TPZAutoPointer<TPZCompMesh> cmesh =
        hcurltest::CreateCMesh(meshType,pOrder);
    SECTION("Vector traces"+MMeshType_Name(meshType)){
        hcurltest::TestVectorTracesUniformMesh(cmesh,meshType);
    }
    for(int k = 1; k < maxK; k++){
        SECTION("Funcion traces "+MMeshType_Name(meshType)+" p"+std::to_string(k)){
            hcurltest::TestFunctionTracesUniformMesh(cmesh,meshType,k);
        }
    }
}

TEST_CASE("Testing curl of HCurl functions",
          "[hcurl_tests][mesh]") {
    constexpr int pOrder{1};
    constexpr int maxK{5};
    auto meshType = GENERATE(MMeshType::ETriangular,
                             MMeshType::EQuadrilateral,
                             MMeshType::ETetrahedral,
                             MMeshType::EHexahedral
                             // MMeshType::EPrismatic//NEEDSFIX
                             );
    TPZAutoPointer<TPZCompMesh> cmesh =
        hcurltest::CreateCMesh(meshType,pOrder);
    
    for(int k = 1; k < maxK; k++){
        SECTION("Funcion curl "+MMeshType_Name(meshType)+" p"+std::to_string(k)){
            hcurltest::TestFunctionCurlUniformMesh(cmesh,meshType,k);
        }
    }
}


namespace hcurltest{

    void TestVectorTracesUniformMesh(TPZAutoPointer<TPZCompMesh> cmesh, MMeshType type){
        const int dim = MMeshType_Dimension(type);
        MElementType elType = [&](){
            switch(type){
            case MMeshType::ETriangular: return ETriangle;
            case MMeshType::EQuadrilateral: return EQuadrilateral;
            case MMeshType::ETetrahedral: return ETetraedro;
            case MMeshType::EHexahedral: return ECube;
            case MMeshType::EPrismatic: return EPrisma;
            case MMeshType::EPyramidal: return EPiramide;
            default: return ENoType;
            }
        }();
        for(auto dummyCel : cmesh->ElementVec()){
            const auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
            const auto gel = cel->Reference();
            //skips boundary els
            if(!cel || cel->Reference()->Type() != elType) continue;

            TPZMaterialDataT<STATE> elData;
            cel->InitMaterialData(elData);
            const int elNNodes = MElementType_NNodes(elType);
            const auto *gmesh =
                cmesh->Reference();
            for (auto iCon = 0; iCon <cel->NConnects(); iCon++) {
                auto &con = cel->Connect(iCon);

                if(con.NElConnected() < 2) continue;

                const int iSide = iCon + elNNodes;

                TPZTransform<> elTransform(gel->SideToSideTransform(iSide, gel->NSides() - 1));
                TPZGeoElSide gelSide(gel, iSide);
                const auto sideDim = gelSide.Dimension();
                //the following vector will be the edge tg vector if 2D, the normal vector if 3D
                TPZManVector<REAL,3> vec(3,0);
                switch(sideDim){
                case 1:{
                    TPZManVector<int, 2> edgeNodes(2, 0);
                    for (auto i = 0; i < 2; i++) edgeNodes[i] = gel->SideNodeIndex(iSide, i);
                    const REAL sign = edgeNodes[0] < edgeNodes[1] ? 1 : -1;
                    REAL edgeLength = 0;
                    TPZVec<REAL> p0(3), p1(3);
                    gmesh->NodeVec()[edgeNodes[0]].GetCoordinates(p0);
                    gmesh->NodeVec()[edgeNodes[1]].GetCoordinates(p1);
                    edgeLength =
                        std::sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                                  (p0[2] - p1[2]) * (p0[2] - p1[2]));
                    for (auto x = 0; x < 3; x++) vec[x] = sign * (p1[x] - p0[x]) / edgeLength;
                }
                    break;
                case 2:{
                    TPZManVector<REAL,3> xCenter(2,0);
                    switch(gelSide.NSideNodes()){
                    case 3:
                        pztopology::TPZTriangle::CenterPoint(pztopology::TPZTriangle::NSides-1,xCenter);
                        break;
                    case 4:
                        pztopology::TPZQuadrilateral::CenterPoint(pztopology::TPZQuadrilateral::NSides-1,xCenter);
                        break;
                    default:
                        DebugStop();
                    }
                    gelSide.Normal(xCenter,vec);
                }
                    break;
                default:
                    DebugStop();
                }

                TPZManVector<REAL,3> centerPt(sideDim,0);
                gelSide.CenterPoint(centerPt);
                TPZManVector<REAL,3> ptEl(gel->Dimension(),0);
                elTransform.Apply(centerPt,ptEl);
                cel->ComputeRequiredData(elData, ptEl);
                TPZFNMatrix<60,REAL> elDeformedDirections(elData.fDeformedDirections);

                //gather all the sides contained in the closure of iSide
                TPZStack<int> smallSides;
                gel->LowerDimensionSides(iSide,smallSides);
                smallSides.Push(iSide);//include the side itself
                const int pOrder = cmesh->GetDefaultOrder();
                TPZGeoElSide neighGelSide = gelSide.Neighbour();
                while(neighGelSide != gelSide) {
                    const auto neighCel = dynamic_cast<TPZInterpolatedElement *> (neighGelSide.Element()->Reference());
                    if (!neighCel) {
                        neighGelSide = neighGelSide.Neighbour();
                        continue;
                    }
                    const auto neighGel = neighCel->Reference();
                    const int neighNNodes = MElementType_NNodes(neighGel->Type());
                    const auto neighSide = neighGelSide.Side();
                    const auto neighDim = neighGelSide.Element()->Dimension();
                    TPZTransform<> neighTransform(neighCel->Reference()->SideToSideTransform(neighSide,
                                                                                             neighGel->NSides() -1));
                    TPZTransform<> localTransf(sideDim);
                    gelSide.SideTransform3(neighGelSide,localTransf);
                    TPZMaterialDataT<STATE> neighData;
                    neighCel->InitMaterialData(neighData);

                    TPZManVector <REAL,3> ptsN(neighGelSide.Dimension(),0),ptNeigh(neighDim);
                    localTransf.Apply(centerPt,ptsN);
                    neighTransform.Apply(ptsN,ptNeigh);
                    neighCel->ComputeRequiredData(neighData, ptNeigh);

                    TPZFNMatrix<60,REAL> neighDeformedDirections(neighData.fDeformedDirections);
                    //since this test was based on the basis functions tests, the following std::set
                    //is a lazy way to avoid multiple testing of the same vectors
                    std::set<int64_t> testedVectors;
                    for(auto subSide : smallSides){
                        TPZGeoElSide gelSubSide(gel, subSide);
                        if(gel->SideDimension(subSide) < 1) continue;
                        const int subConnect = subSide - elNNodes;
                        const int neighSubSide = [&](){
                            TPZGeoElSide neighGelSubSide = gelSubSide.Neighbour();
                            while(neighGelSubSide.Element() != neighGelSide.Element()) {
                                neighGelSubSide = neighGelSubSide.Neighbour();
                                if(neighGelSubSide.Element() == gel){
                                    DebugStop();
                                }
                            }
                            return neighGelSubSide.Side();
                        }();
                        const int nShapes = cel->NConnectShapeF(subConnect,pOrder);

                        const int firstElShape = [&](){
                            int firstElShapeTemp = 0;
                            for(auto jCon = 0; jCon < subConnect; jCon++){
                                firstElShapeTemp += cel->NConnectShapeF(jCon,cel->EffectiveSideOrder(jCon+elNNodes));
                            }
                            return firstElShapeTemp;
                        }();
                    

                        const int firstNeighShape = [&](){
                            int firstNeighShapeTemp = 0;
                            const int neighCon = neighSubSide - neighGelSide.Element()->NNodes();
                            for(auto jCon = 0; jCon < neighCon; jCon++){
                                firstNeighShapeTemp += neighCel->NConnectShapeF(jCon,neighCel->EffectiveSideOrder(jCon+neighNNodes));
                            }
                            return firstNeighShapeTemp;
                        }();

                        TPZManVector<REAL,3> elVec(3,0), neighVec(3,0);
                        for(auto iShape = 0; iShape < nShapes; iShape ++){
                            const int elPhiIndex = firstElShape+iShape;
                            const int neighPhiIndex = firstNeighShape+iShape;
                            const int elVecIndex = elData.fVecShapeIndex[elPhiIndex].first;
                            const int neighVecIndex = neighData.fVecShapeIndex[neighPhiIndex].first;
                            if (testedVectors.find(elVecIndex) == testedVectors.end()) {
                                testedVectors.insert(elVecIndex);
                            }else{
                                continue;
                            }
                            for (auto x = 0; x < 3; x++) elVec[x] = elDeformedDirections(x,elVecIndex);
                            for (auto x = 0; x < 3; x++) neighVec[x] = neighDeformedDirections(x,neighVecIndex);
                            REAL diffTrace{0};
                            TPZManVector<REAL,3> elTrace,neighTrace;
                            const bool checkTraces = CheckTracesFunc(diffTrace,elVec,neighVec,vec,sideDim,elTrace,neighTrace);
                            CAPTURE(gel->Index(),neighGel->Index());
                            CAPTURE(elVecIndex,neighVecIndex);
                            CAPTURE(elData.jacobian,neighData.jacobian);
                            REQUIRE(checkTraces);
                        }
                    }
                    neighGelSide = neighGelSide.Neighbour();
                }
            }
        }
    }//hcurltest::TestVectorTracesUniformMesh

    void TestFunctionTracesUniformMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                                       MMeshType type, const int pOrder){
        const int dim = MMeshType_Dimension(type);   
        MElementType elType = [&](){
            switch(type){
            case MMeshType::ETriangular: return ETriangle;
            case MMeshType::EQuadrilateral: return EQuadrilateral;
            case MMeshType::ETetrahedral: return ETetraedro;
            case MMeshType::EHexahedral: return ECube;
            case MMeshType::EPrismatic: return EPrisma;
            case MMeshType::EPyramidal: return EPiramide;
            default: return ENoType;
            }
        }();

        if(pOrder > 1 ){
            cmesh->SetDefaultOrder(pOrder);
            for(auto cel : cmesh->ElementVec()){
                TPZInterpolatedElement *intel =
                    dynamic_cast<TPZInterpolatedElement *> (cel);
                if(intel) intel->PRefine(pOrder);
            }
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();
        }
        const auto *gmesh =
            cmesh->Reference();
        for(auto dummyCel : cmesh->ElementVec()){
            const auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
            if(!cel) continue;
            const auto gel = cel->Reference();
            //skips boundary els
            if(!cel || cel->Reference()->Type() != elType) continue;

            TPZMaterialDataT<STATE> elData;
            cel->InitMaterialData(elData);
            const int nState = cel->Material()->NStateVariables();
            const int elNNodes = MElementType_NNodes(elType);
            for (auto iCon = 0; iCon <cel->NConnects(); iCon++) {
                auto &con = cel->Connect(iCon);
                const int nShape = con.NShape();
                {
                    const bool check = con.NDof( *cmesh) == nShape * nState;
                    REQUIRE(check);
                }

                if(con.NElConnected() < 2) continue;

                const int iSide = iCon + elNNodes;
                const int pOrderIntRule = cel->EffectiveSideOrder(iSide)*2;
                TPZIntPoints *sideIntRule = gel->CreateSideIntegrationRule(iSide, pOrderIntRule);
                const int npts = sideIntRule->NPoints();

                TPZTransform<> elTransform(gel->SideToSideTransform(iSide, gel->NSides() - 1));
                TPZGeoElSide gelSide(gel, iSide);
                const auto sideDim = gelSide.Dimension();

                //the following vector will be the edge tg vector if 2D, the normal vector if 3D
                TPZManVector<REAL,3> vec(3,0);
                switch(sideDim){
                case 1:{
                    TPZManVector<int, 2> edgeNodes(2, 0);
                    for (auto i = 0; i < 2; i++) edgeNodes[i] = gel->SideNodeIndex(iSide, i);
                    const REAL sign = edgeNodes[0] < edgeNodes[1] ? 1 : -1;
                    REAL edgeLength = 0;
                    TPZVec<REAL> p0(3), p1(3);
                    gmesh->NodeVec()[edgeNodes[0]].GetCoordinates(p0);
                    gmesh->NodeVec()[edgeNodes[1]].GetCoordinates(p1);
                    edgeLength =
                        std::sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                                  (p0[2] - p1[2]) * (p0[2] - p1[2]));
                    for (auto x = 0; x < 3; x++) vec[x] = sign * (p1[x] - p0[x]) / edgeLength;
                }
                    break;
                case 2:{
                    TPZManVector<REAL,3> xCenter(2,0);
                    switch(gelSide.NSideNodes()){
                    case 3:
                        pztopology::TPZTriangle::CenterPoint(pztopology::TPZTriangle::NSides-1,xCenter);
                        break;
                    case 4:
                        pztopology::TPZQuadrilateral::CenterPoint(pztopology::TPZQuadrilateral::NSides-1,xCenter);
                        break;
                    default:
                        DebugStop();
                    }
                    gelSide.Normal(xCenter,vec);
                }
                    break;
                default:
                    DebugStop();
                }
                bool firstNeighbour{true};


                //gather all the sides contained in the closure of iSide
                TPZStack<int> smallSides;
                gel->LowerDimensionSides(iSide,smallSides);
                smallSides.Push(iSide);//include the side itself
                //iterate on the neighbouring sides
                TPZGeoElSide neighGelSide = gelSide.Neighbour();
                while(neighGelSide != gelSide) {
                    const auto neighCel = dynamic_cast<TPZInterpolatedElement *> (neighGelSide.Element()->Reference());
                    if (!neighCel) {
                        neighGelSide = neighGelSide.Neighbour();
                        continue;
                    }
                    const auto neighGel = neighCel->Reference();
                    const int neighNNodes = MElementType_NNodes(neighGel->Type());
                    const auto neighSide = neighGelSide.Side();
                    TPZTransform<> neighTransform(neighCel->Reference()->SideToSideTransform(neighSide,
                                                                                             neighCel->Reference()->NSides() -1));
                    TPZTransform<> localTransf(sideDim);
                    gelSide.SideTransform3(neighGelSide,localTransf);
                    TPZMaterialDataT<STATE> neighData;
                    neighCel->InitMaterialData(neighData);

                    const auto neighDim = neighGelSide.Element()->Dimension();
                    TPZManVector <REAL,3> pts(sideDim,0),ptsN(sideDim,0), ptEl(dim,0),ptNeigh(neighDim,0);
                    REAL w;
                    TPZFNMatrix<60,REAL> elShape,neighShape;

                    for(auto subSide : smallSides){
                        TPZGeoElSide gelSubSide(gel, subSide);
                        if(gel->SideDimension(subSide) < 1) continue;
                        const int subConnect = subSide - elNNodes;
                        const int neighSubSide = [&](){
                            TPZGeoElSide neighGelSubSide = gelSubSide.Neighbour();
                            while(neighGelSubSide.Element() != neighGelSide.Element()) {
                                neighGelSubSide = neighGelSubSide.Neighbour();
                                if(neighGelSubSide.Element() == gel){
                                    DebugStop();
                                }
                            }
                            return neighGelSubSide.Side();
                        }();

                        const int firstElShape = [&](){
                            int firstElShapeTemp = 0;
                            for(auto jCon = 0; jCon < subConnect; jCon++){
                                firstElShapeTemp += cel->NConnectShapeF(jCon,cel->EffectiveSideOrder(jCon+elNNodes));
                            }
                            return firstElShapeTemp;
                        }();
                        const int nShapes = cel->NConnectShapeF(subConnect,pOrder);

                        {
                            const auto elConIndex = cel->SideConnectIndex(cel->NSideConnects(subSide) - 1, subSide);
                            const auto neighConIndex = neighCel->SideConnectIndex(neighCel->NSideConnects(neighSubSide) - 1, neighSubSide);
                            const bool check = elConIndex == neighConIndex;
                            REQUIRE(check);
                        }
                        const int firstNeighShape = [&](){
                            int firstNeighShapeTemp = 0;
                            const int neighCon = neighSubSide - neighGelSide.Element()->NNodes();
                            for(auto jCon = 0; jCon < neighCon; jCon++){
                                firstNeighShapeTemp += neighCel->NConnectShapeF(jCon,neighCel->EffectiveSideOrder(jCon+neighNNodes));
                            }
                            return firstNeighShapeTemp;
                        }();

                        for (auto ipt = 0; ipt < npts; ipt++) {
                            sideIntRule->Point(ipt, pts, w);

                            elTransform.Apply(pts,ptEl);
                            cel->ComputeRequiredData(elData, ptEl);
                            TPZHCurlAuxClass::ComputeShape(elData.fVecShapeIndex, elData.phi,
                                                           elData.fDeformedDirections,elShape);
                            localTransf.Apply(pts,ptsN);
                            neighTransform.Apply(ptsN,ptNeigh);
                            neighCel->ComputeRequiredData(neighData, ptNeigh);
                            TPZHCurlAuxClass::ComputeShape(neighData.fVecShapeIndex, neighData.phi,
                                                           neighData.fDeformedDirections,neighShape);

                            TPZManVector<REAL,3> elShapeFunc(3,0), neighShapeFunc(3,0);
                            bool anyWrongCheck = false;
                            for(auto iShape = 0; iShape < nShapes; iShape ++){
                                const int elPhiIndex = firstElShape+iShape;
                                const int neighPhiIndex = firstNeighShape+iShape;
                                const int elH1phiIndex = elData.fVecShapeIndex[elPhiIndex].second;
                                const int neighH1phiIndex = neighData.fVecShapeIndex[neighPhiIndex].second;
                                const bool checkPhis = std::abs(elData.phi(elH1phiIndex,0) - neighData.phi(neighH1phiIndex,0)) < tol;

                                REQUIRE(checkPhis);
                                anyWrongCheck = !checkPhis || anyWrongCheck;
                                if(anyWrongCheck) {
                                    break;
                                }
                                for (auto x = 0; x < 3; x++) elShapeFunc[x] = elShape(elPhiIndex,x);
                                for (auto x = 0; x < 3; x++) neighShapeFunc[x] = neighShape(neighPhiIndex,x);
                                REAL diffTrace{0};
                                TPZManVector<REAL,3> elTrace,neighTrace;
                                const bool checkTraces = CheckTracesFunc(diffTrace,elShapeFunc,neighShapeFunc,vec,sideDim,elTrace,neighTrace);
                                
                                CAPTURE(gel->Index(),neighGel->Index());                                
                                CAPTURE(iSide,neighSide);
                                CAPTURE(subSide,neighSubSide);
                                CAPTURE(ptEl,elData.x);
                                CAPTURE(ptNeigh,neighData.x);
                                CAPTURE(elPhiIndex,neighPhiIndex);
                                CAPTURE(elShapeFunc,neighShapeFunc);
                                REQUIRE(checkTraces);

                                anyWrongCheck = !checkTraces || anyWrongCheck;
                                if(anyWrongCheck) {
                                    break;
                                }
                            }

                            //now checking sideshapefunction. this needs to be checked only once
                            if(firstNeighbour){
                                TPZFNMatrix<30,REAL> sideShapeFuncs,sideShapeCurl;
                                cel->SideShapeFunction(gelSide.Side(),pts,sideShapeFuncs,sideShapeCurl);

                                const int firstSideShape = [&](){
                                    int firstSideShapeTemp = 0;
                                    TPZStack<int> subSubSides;
                                    gel->LowerDimensionSides(gelSide.Side(),subSubSides);
                                    for(auto subsubSide : subSubSides){
                                        if(gel->SideDimension(subsubSide) < 1) continue;
                                        if(subsubSide == subSide) break;
                                        firstSideShapeTemp += cel->NConnectShapeF(subsubSide - elNNodes,pOrder);
                                    }
                                    return firstSideShapeTemp;
                                }();
                                TPZManVector<REAL,3> sideShapeFunc(3,0);
                                for(auto iShape = 0; iShape < nShapes; iShape ++){
                                    for (auto x = 0; x < 3; x++) elShapeFunc[x] = elShape(firstElShape + iShape,x);
                                    for (auto x = 0; x < 3; x++) sideShapeFunc[x] = sideShapeFuncs(firstSideShape + iShape,x);
                                    REAL diffTrace{0};
                                    TPZManVector<REAL,3> elTrace, sideTrace;
                                    const bool checkTraces = CheckTracesFunc(diffTrace,elShapeFunc,sideShapeFunc,vec,sideDim,elTrace,sideTrace);
                                    CAPTURE(elShapeFunc,sideShapeFunc);
                                    CAPTURE(elTrace,sideTrace);
                                    REQUIRE(checkTraces);
                                    anyWrongCheck = !checkTraces || anyWrongCheck;
                                }
                                if(anyWrongCheck){
                                    break;
                                }
                            }
                        }
                    }
                    firstNeighbour = false;
                    neighGelSide = neighGelSide.Neighbour();
                }
            }
        }
        /***********************************************************************************************************
         *              the following lines might be useful for analysing the basis functions
         ***********************************************************************************************************/

        // TPZLinearAnalysis an(cmesh, false);
        // const int postProcessResolution = 3;
        // const std::string executionInfo = [&]() {
        //   std::string name("");
        //   name.append(MMeshType_Name(type));
        //   name.append(std::to_string(pOrder));
        //   return name;
        // }();

        // const std::string plotfile =
        //     "solution" + executionInfo + ".vtk"; // where to print the vtk files
        // TPZStack<std::string> scalnames, vecnames;
        // vecnames.Push("Solution"); // print the state variable
        // TPZFMatrix<STATE> sol(an.Solution().Rows(), 1);
        // sol.Zero();
        // for (int i = 0; i < sol.Rows(); i++) {
        //   sol(i - 1 < 0 ? 0 : i - 1, 0) = 0;
        //   sol(i, 0) = 1;
        //   an.LoadSolution(sol);
        //   an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        //   an.SetStep(i);
        //   an.PostProcess(postProcessResolution);
        // }
    }//hcurltest::TestFunctionTracesUniformMesh

    void TestFunctionCurlUniformMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                                       MMeshType type, const int pOrder)
    {
        const auto oldPrecision = Catch::StringMaker<REAL>::precision;
        Catch::StringMaker<REAL>::precision = std::numeric_limits<REAL>::max_digits10;
        const int dim = MMeshType_Dimension(type);
        const int curlDim = dim == 1 ? 1 : 2 * dim -3;
        MElementType elType = [&](){
            switch(type){
            case MMeshType::ETriangular: return ETriangle;
            case MMeshType::EQuadrilateral: return EQuadrilateral;
            case MMeshType::ETetrahedral: return ETetraedro;
            case MMeshType::EHexahedral: return ECube;
            case MMeshType::EPrismatic: return EPrisma;
            case MMeshType::EPyramidal: return EPiramide;
            default: return ENoType;
            }
        }();

        if(pOrder > 1 ){
            cmesh->SetDefaultOrder(pOrder);
            for(auto cel : cmesh->ElementVec()){
                TPZInterpolatedElement *intel =
                    dynamic_cast<TPZInterpolatedElement *> (cel);
                if(intel) intel->PRefine(pOrder);
            }
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();
        }
        
        for(auto dummyCel : cmesh->ElementVec()){
            const auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
            if(!cel) continue;
            const auto gel = cel->Reference();
            //skips boundary els
            if(!cel || cel->Reference()->Type() != elType) continue;

            TPZMaterialDataT<STATE> elData;
            cel->InitMaterialData(elData);
            const int nState = cel->Material()->NStateVariables();
            const int nShapeFunc = cel->NShapeF();
            const int nSides = gel->NSides();
            TPZIntPoints *intRule =
                gel->CreateSideIntegrationRule(nSides-1, pOrder);
            const int nIntPts = intRule->NPoints();
            TPZManVector<REAL,3> qsi(dim);
            REAL w;

            
            auto Cross = [dim,curlDim](const TPZFMatrix<REAL> v1,
                               const TPZFMatrix<REAL>&v2,
                               TPZFMatrix<REAL> &res)
            {
                res.Redim(curlDim,1);
                switch(dim){
                case 2:
                    res(0,0) = v1(0,0) * v2(1,0) - v1(1,0) * v2(0,0);
                    break;
                case 3:
                    res(0,0) = v1(1,0) * v2(2,0) - v1(2,0) * v2(1,0);
                    res(1,0) = v1(2,0) * v2(0,0) - v1(0,0) * v2(2,0);
                    res(2,0) = v1(0,0) * v2(1,0) - v1(1,0) * v2(0,0);
                    break;
                }
            };
            
            TPZFNMatrix<3,REAL>
                dir(dim,1,0),
                gphiHat(dim,1,0),              
                curlX(curlDim,1,0),
                curlCalc(curlDim,1,0);

            
            const auto tol = std::numeric_limits<REAL>::epsilon()*10;
            const auto &indices = elData.fVecShapeIndex;
            const auto &directions = elData.fMasterDirections;
            const auto &axes = elData.axes;
            
            for(auto iPt = 0; iPt < nIntPts; iPt++){
                intRule->Point(iPt, qsi,w);
                cel->ComputeRequiredData(elData, qsi);
                const auto &gradPhiVec = elData.dphi;
                const auto &jac = elData.jacobian;
                const auto &detjac = elData.detjac;
                
                const auto &curlCalcVec = elData.curlphi;
                for(auto iShape = 0; iShape < nShapeFunc; iShape++){
                    curlCalcVec.GetSub(0, iShape, curlDim, 1, curlCalc);
                    const auto &vi = indices[iShape].first;
                    const auto &pi = indices[iShape].second;
                    directions.GetSub(0, vi, dim, 1, dir);
                    gradPhiVec.GetSub(0, pi, dim, 1, gphiHat);
                    Cross(gphiHat, dir ,curlX);
                    if (dim == 3){
                        curlX = jac * curlX;
                    }
                    curlX *= 1./detjac;
                    CAPTURE(iShape, jac, gphiHat, dir, curlX,curlCalc);
                    const auto diff = fabs(Norm(curlX-curlCalc));
                    REQUIRE(diff == Approx(0.0).margin(tol));
                }//for iShape
            }//for iPt
        }//for dummyCel
        Catch::StringMaker<REAL>::precision = oldPrecision;
    }//hcurltest::TestFunctionCurlUniformMesh


    template <class TGEOM>
    void PrintShapeFunctions(const int pOrder){
        constexpr auto dim = TGEOM::Dimension;
        TPZGeoMesh *gmesh =
            TPZGeoMeshTools::CreateGeoMeshSingleElT<TGEOM>(1,false);
        TPZCompMesh *cmesh = new TPZCompMesh (gmesh);

        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(dim);

        auto mat = new TPZHCurlProjection(1,dim);
        cmesh->InsertMaterialObject(mat);
        cmesh->SetAllCreateFunctionsHCurl();
        cmesh->AutoBuild();
        cmesh->CleanUpUnconnectedNodes();

        TPZLinearAnalysis an(cmesh,false);
        const int postProcessResolution = 3;
        const std::string executionInfo = [&](){
            std::string name("");
            name.append(MElementType_Name(TGEOM::Type()));
            name.append(std::to_string(pOrder));
            return name;
        }();

        const std::string plotfile = "shapeFuncs"+executionInfo+".vtk";//where to print the vtk files
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("Solution");//print the state variable
        TPZFMatrix<STATE>& sol = an.Solution();
        sol.Zero();
        for(int i = 0; i < sol.Rows(); i++){
            sol(i - 1 < 0 ? 0 : i - 1 , 0) = 0;
            sol(i,0) = 1;
            an.LoadSolution(sol);
            an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
            an.SetStep(i);
            an.PostProcess(postProcessResolution);
        }
        delete cmesh;
        delete gmesh;
    }


    TPZAutoPointer<TPZCompMesh> CreateCMesh(MMeshType type, int pOrder){
        constexpr int ndiv{1};
        const int dim = MMeshType_Dimension(type);
        TPZManVector<int,3> nDivVec(dim,ndiv);
        TPZManVector<REAL,3> minX({0,0,0});
        TPZManVector<REAL,3> maxX({1,1,1});
        if(dim == 2) maxX[2] = 0.0;
            
        TPZManVector<int,2> matIds(2*dim+1,-1);
        matIds[0]=1;
        constexpr bool createBoundEls{true};
        TPZAutoPointer<TPZGeoMesh>gmesh =
            TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX,maxX,matIds,nDivVec,type,createBoundEls);
        //            {
        //                TPZVec<TPZGeoEl*> sons;
        //                gmesh->Element(0)->Divide(sons);
        //            }
#ifdef HCURL_OUTPUT_TXT
        std::ofstream outTXT("gmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
        gmesh->Print(outTXT);
        outTXT.close();
#endif
#ifdef HCURL_OUTPUT_VTK
        std::ofstream outVTK("gmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK);
        outVTK.close();
#endif
        TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(dim);

        auto mat = new TPZHCurlProjection(matIds[0],dim);
        cmesh->InsertMaterialObject(mat);
        auto bcond = mat->CreateBC(mat,matIds[1],0,TPZFNMatrix<1,STATE>(1,1,0),{0});
        cmesh->InsertMaterialObject(bcond);
        cmesh->SetAllCreateFunctionsHCurl();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
    }
            
    bool CheckTracesFunc(REAL &diffTrace, const TPZVec<REAL> &elShapeFunc, const TPZVec<REAL>&neighShapeFunc,
                         const TPZVec<REAL> &vec, const int sideDim, TPZVec<REAL> &elTrace, TPZVec<REAL> &neighTrace){
        REAL averageTrace{0};
        diffTrace = 0;
        switch(sideDim){
        case 1:{
            elTrace.Resize(1,0);
            neighTrace.Resize(1,0);
            for (auto x = 0; x < 3; x++) elTrace[0] += elShapeFunc[x]*vec[x];
            for (auto x = 0; x < 3; x++) neighTrace[0] += neighShapeFunc[x]*vec[x];
            diffTrace = std::abs(elTrace[0] - neighTrace[0]);
            averageTrace = (sqrt(elTrace[0]*elTrace[0])+sqrt(neighTrace[0]*neighTrace[0]))/2;
        }
            break;
        case 2:{
            TPZManVector<REAL,3> temp(3,0);
            elTrace.Resize(3,0);
            neighTrace.Resize(3,0);
            VectorProduct(vec,elShapeFunc,temp);
            VectorProduct(temp,vec,elTrace);

            VectorProduct(vec,neighShapeFunc,temp);
            VectorProduct(temp,vec,neighTrace);
            for(auto x = 0; x < 3; x++) diffTrace+= (neighTrace[x]-elTrace[x])*(neighTrace[x]-elTrace[x]);
            diffTrace= sqrt(diffTrace);
            REAL firstTrace = 0, secondTrace = 0;
            for(auto x = 0; x < 3; x++) firstTrace+= elTrace[x]*elTrace[x];
            for(auto x = 0; x < 3; x++) secondTrace+= neighTrace[x]*neighTrace[x];
            firstTrace = sqrt(firstTrace);
            secondTrace = sqrt(secondTrace);
            averageTrace =(firstTrace+secondTrace)/2;
        }
            break;
        default:
            return false;
        }
        if(averageTrace < tol)  averageTrace = 1;

        return  diffTrace/averageTrace < tol;
    }

    void ComputeDirections (TPZGeoEl *gel, TPZFMatrix<REAL> &deformedDirections,
                            TPZFMatrix<REAL> &masterDirections, const TPZVec<int> &transformationIds) {
        int dim = gel->Dimension();
        TPZFMatrix<REAL> gradX(dim, dim, 0);
        for (auto x = 0; x < dim; x++) gradX(x, x) = 1;
        const int nVec = dim * gel->NSides();
        masterDirections.Resize(dim, nVec);
        deformedDirections.Resize(3, nVec);
        switch (gel->Type()) {
        case EOned:
            pztopology::TPZLine::ComputeHCurlDirections(gradX,
                                                        masterDirections,transformationIds);//these are the vectors on the master element
            break;
        case ETriangle:
            pztopology::TPZTriangle::ComputeHCurlDirections(gradX,
                                                            masterDirections,transformationIds);//these are the vectors on the master element
            break;
        case EQuadrilateral:
            pztopology::TPZQuadrilateral::ComputeHCurlDirections(gradX,
                                                                 masterDirections,transformationIds);//these are the vectors on the master element
            break;
        case ETetraedro:
            pztopology::TPZTetrahedron::ComputeHCurlDirections(gradX,
                                                               masterDirections,transformationIds);//these are the vectors on the master element
            break;
        case EPiramide:
            DebugStop();//pztopology::TPZPyramid::ComputeHCurlDirections(gradX, masterDirections,transformationIds);//these are the vectors on the master element
            break;
        case EPrisma:
            pztopology::TPZPrism::ComputeHCurlDirections(gradX, masterDirections,transformationIds);//these are the vectors on the master element
            break;
        case ECube:
            pztopology::TPZCube::ComputeHCurlDirections(gradX,
                                                        masterDirections,transformationIds);//these are the vectors on the master element
            break;
        default:
            DebugStop();
        }
        TPZVec<REAL> xiCenter(dim, 0);
        TPZFMatrix<REAL> jacobian(dim, dim, 0), axes(dim, 3), jacInv(dim, dim, 0);
        REAL detJac = 0;
        gel->Jacobian(xiCenter, jacobian, axes, detJac, jacInv);

#ifdef VERBOSE_HCURL
        std::cout << '\n';
        jacobian.Print("Original Jacobian:");
        axes.Print("Original Axes:");
#endif
        for (auto iVec = 0; iVec < nVec; iVec++) {

            TPZManVector<REAL, 3> tempDirection(dim, 0);
            for (auto i = 0; i < dim; i++) {
                //covariant piola transform: J^{-T}
                tempDirection[i] = 0;
                for (auto j = 0; j < dim; j++) tempDirection[i] += jacInv(j, i) * masterDirections(j, iVec);
            }
            for (auto i = 0; i < 3; i++) {
                deformedDirections(i, iVec) = 0;
                for (auto j = 0; j < dim; j++) deformedDirections(i, iVec) += axes(j, i) * tempDirection[j];
            }
        }
    }//ComputeDirections

    void VectorProduct(const TPZVec<REAL> &v1, const TPZVec<REAL> &v2,TPZVec<REAL> &result){
        REAL x1=v1[0], y1=v1[1],z1=v1[2];
        REAL x2=v2[0], y2=v2[1],z2=v2[2];
        result.Resize(v1.NElements());
        result[0]=y1*z2-z1*y2;
        result[1]=z1*x2-x1*z2;
        result[2]=x1*y2-y1*x2;
    };//VectorProduct
}//namespace

#ifdef VERBOSE_HCURL
#undef VERBOSE_HCURL
#endif
#ifdef HCURL_OUTPUT_TXT
#undef HCURL_OUTPUT_TXT
#endif
#ifdef HCURL_OUTPUT_VTK
#undef HCURL_OUTPUT_VTK
#endif