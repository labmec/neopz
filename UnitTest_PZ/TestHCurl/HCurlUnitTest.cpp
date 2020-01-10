/**
 * @file HCurlUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of HCurl-conforming approximation spaces
 *
 */
#include <iostream>

#include "pzlog.h"
#include "TPZTopologyUtils.h"
#include "pzgmesh.h"

#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

#include "pzgeoelrefless.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testhcurl"));
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

#include "pzgengrid.h"
#include "TPZGenGrid3D.h"
#include "pzcmesh.h"
#include "TPZMatHelmholtz.cpp"
#include "pzanalysis.h"
#include "pzintel.h"
#include "TPZCompElHCurl.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZVTKGeoMesh.h"
// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz #define BOOST_TEST_MAIN pz hcurl_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>


#define NOISY_HCURL //outputs useful debug info
//#define NOISY_HCURL_VTK


struct SuiteInitializer{
    SuiteInitializer(){
        InitializePZLOG();
        boost::unit_test::unit_test_log.set_threshold_level( boost::unit_test::log_warnings );
    }
};
BOOST_FIXTURE_TEST_SUITE(hcurl_tests,SuiteInitializer)
    namespace hcurltest{
        constexpr REAL tol = 1e-10;
        template <class TTopol>
        void CompareVectorTraces(const TPZFMatrix<REAL> &);
        
        void TestTracesUniformMesh(MElementType type, const int pOrder, const int dim);

        template <class TGEOM>
        void PrintShapeFunctions(const int pOrder);

        namespace auxiliaryfuncs{
            void ComputeDirections (TPZGeoEl *, TPZFMatrix<REAL> &, TPZFMatrix<REAL> &,const TPZVec<int> &);

            void VectorProduct(const TPZVec<REAL> &, const TPZVec<REAL> &, TPZVec<REAL> &);

            void ComputeNormal(TPZGeoMesh *, TPZVec<int64_t>, TPZVec<REAL> &);

            TPZGeoMesh *CreateGeoMesh2D(int nelx, int nely, MElementType meshType, TPZVec<int> &matIds);

            TPZGeoMesh *CreateGeoMesh3D(int nelx, int nely, int nelz, MElementType meshType, TPZVec<int> &matIds);
        }
    }


    BOOST_AUTO_TEST_CASE(hcurl_topology_tests) {
        {
            TPZFMatrix<REAL> nodeCoords(3,3);
            nodeCoords(0,0) = -1;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  0;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            hcurltest::CompareVectorTraces<pztopology::TPZTriangle>(nodeCoords);
        }
        {
            TPZFMatrix<REAL> nodeCoords(4,3);
            nodeCoords(0,0) =  0;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  1;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            nodeCoords(3,0) =  0;   nodeCoords(3,1) =  1;   nodeCoords(3,2) =  0;
            hcurltest::CompareVectorTraces<pztopology::TPZQuadrilateral>(nodeCoords);
        }

        {
            TPZFMatrix<REAL> nodeCoords(4,3);
            nodeCoords(0,0) = -1;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  0;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            nodeCoords(3,0) =  0;   nodeCoords(3,1) =  0;   nodeCoords(3,2) =  1;
            hcurltest::CompareVectorTraces<pztopology::TPZTetrahedron>(nodeCoords);
        }

        {
            TPZFMatrix<REAL> nodeCoords(8,3);
            nodeCoords(0,0) =  0;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  1;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            nodeCoords(3,0) =  0;   nodeCoords(3,1) =  1;   nodeCoords(3,2) =  0;
            nodeCoords(4,0) =  0;   nodeCoords(4,1) =  0;   nodeCoords(4,2) =  1;
            nodeCoords(5,0) =  1;   nodeCoords(5,1) =  0;   nodeCoords(5,2) =  1;
            nodeCoords(6,0) =  1;   nodeCoords(6,1) =  1;   nodeCoords(6,2) =  1;
            nodeCoords(7,0) =  0;   nodeCoords(7,1) =  1;   nodeCoords(7,2) =  1;
            hcurltest::CompareVectorTraces<pztopology::TPZCube>(nodeCoords);
        }

        {
            TPZFMatrix<REAL> nodeCoords(6,3);
            nodeCoords(0,0) = -1;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  0;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            nodeCoords(3,0) = -1;   nodeCoords(3,1) =  0;   nodeCoords(3,2) =  1;
            nodeCoords(4,0) =  1;   nodeCoords(4,1) =  0;   nodeCoords(4,2) =  1;
            nodeCoords(5,0) =  0;   nodeCoords(5,1) =  1;   nodeCoords(5,2) =  1;
            hcurltest::CompareVectorTraces<pztopology::TPZPrism>(nodeCoords);
        }
    }

    BOOST_AUTO_TEST_CASE(hcurl_mesh_tests) {
//        hcurltest::PrintShapeFunctions<pzgeom::TPZGeoTetrahedra>(4);
        constexpr int dim2D{2},dim3D{3};
        for(auto k = 1; k <= 5; k++) hcurltest::TestTracesUniformMesh(ETriangle, k,dim2D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestTracesUniformMesh(EQuadrilateral, k,dim2D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestTracesUniformMesh(ETetraedro, k,dim3D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestTracesUniformMesh(ECube, k,dim3D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestTracesUniformMesh(EPrisma, k,dim3D);

    }

    namespace hcurltest{


        template <class TTopol>
        void CompareVectorTraces(const TPZFMatrix<REAL> &nodeCoords) {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            auto elType = TTopol::Type();

            const auto nSides = TTopol::NumSides();
            const auto nNodes = TTopol::NumSides(0);
            const auto nEdges = TTopol::NumSides(1);
            const auto nFaces = TTopol::NumSides(2);
            constexpr auto dim = TTopol::Dimension;
            //mesh
            const auto nPermutations = TTopol::NPermutations;
            TPZManVector<int,21> permVec(nSides,0);
            for(auto iPerm = 0; iPerm< nPermutations; iPerm ++){
                std::cout<<"\tPermutation "<<iPerm+1<<" out of "<<nPermutations<<std::endl;
                pztopology::GetPermutation<TTopol>(iPerm,permVec);
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                for (int iNode = 0; iNode < nNodes; iNode++) {
                    TPZManVector<REAL,21> xNode(3,0);
                    for (auto x = 0; x < 3; x++) xNode[x] = nodeCoords.GetVal(permVec[iNode], x);
                    auto newIndex = gmesh->NodeVec().AllocateNewElement();
                    gmesh->NodeVec()[newIndex].Initialize(xNode, *gmesh);
                }

                TPZVec<int64_t> nodes(nNodes, 0);
                for (auto i = 0; i < nNodes; i++) nodes[i] = permVec[i];
                int64_t index;
                TPZGeoEl *gel = gmesh->CreateGeoElement(elType, nodes, 1, index, 0);
                //computing transformation id for sides
                TPZManVector<int,21> transformationIds(nSides - nNodes,-1);
                for(auto iSide = 0 ; iSide < nEdges+nFaces + dim - 2; iSide++){
                    transformationIds[iSide] = TTopol::GetTransformId(nNodes + iSide, nodes);
                }
                //computing directions
                TPZFMatrix<REAL> masterDirections(0,0), deformedDirections(0,0);
                auxiliaryfuncs::ComputeDirections(gel,deformedDirections,masterDirections,transformationIds);
                //calculating tangent vectors for all edges
                TPZFMatrix<REAL> edgeTangentVecs(nEdges,3);
                TPZVec<REAL> edgeLengthVec(nEdges,0);
                for (auto iEdge = 0; iEdge < nEdges; iEdge++) {
                    const int edgeIndex = nNodes + iEdge;

                    TPZManVector<int, 2> edgeNodes(2, 0);
                    for (auto i = 0; i < 2; i++) edgeNodes[i] = gel->SideNodeIndex(edgeIndex, i);
                    const REAL sign = edgeNodes[0] < edgeNodes[1] ? 1 : -1;
                    TPZVec<REAL> edgeTgVector(3, 0);
                    REAL edgeLength = 0;
                    {
                        TPZVec<REAL> p0(3), p1(3);
                        gmesh->NodeVec()[edgeNodes[0]].GetCoordinates(p0);
                        gmesh->NodeVec()[edgeNodes[1]].GetCoordinates(p1);
                        edgeLength =
                                std::sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                                          (p0[2] - p1[2]) * (p0[2] - p1[2]));
                        for (auto x = 0; x < 3; x++) edgeTangentVecs(iEdge,x) = sign * (p1[x] - p0[x]) / edgeLength;
                        edgeLengthVec[iEdge] = edgeLength;
                    }
                }
                //testing directions associated with edges
                for (auto iEdge = 0; iEdge < nEdges; iEdge++) {
                    const int edgeIndex = nNodes + iEdge;


                    TPZVec<REAL> edgeTgVector(3, 0);
                    for (auto i = 0; i < 3; i++) edgeTgVector[i] = edgeTangentVecs(iEdge,i);
                    REAL edgeLength = edgeLengthVec[iEdge];


                    REAL trace = 0;
                    for (auto x = 0; x < 3; x++) trace += deformedDirections(x, 2 * nEdges + iEdge) * edgeTgVector[x];
                    trace *= edgeLength;
#ifdef NOISY_HCURL
                    std::cout << "\tedgeIndex " << edgeIndex <<"\ttesting ve vectors" << std::endl;
                    std::cout << "\t\ttangent vector: " << std::endl;
                    for (auto x = 0; x < 3; x++) std::cout << "\t\t" << edgeTgVector[x];
                    std::cout << std::endl;
                    std::cout << "\t\tve vector:" << std::endl;
                    for (auto x = 0; x < 3; x++) std::cout << "\t\t" << deformedDirections(x, 2 * nEdges + iEdge) << "\t";
                    std::cout << std::endl;
                    std::cout << "\t\ttrace:" << trace << std::endl;
                    std::cout << std::endl;
#endif
                    bool testTrace = std::abs(trace - 1) < tol;
                    BOOST_CHECK(testTrace);
                    for(auto iVec = 0; iVec <  2; iVec++) {// 2 v^{e,a}
                        trace = 0;
                        for(auto x = 0; x < 3; x++) trace += edgeTgVector[x] * deformedDirections(x,2*iEdge+iVec);
                        trace *=edgeLength;
#ifdef NOISY_HCURL
                        std::cout<<"\t\tvector  "<<iVec<<" out of "<<2<<std::endl;
                        std::cout<<"\t\t\tdirection "<<iVec<<std::endl;
                        std::cout<<"\t\t\tintegral of tangential trace over edge: "<<trace<<std::endl;
#endif
                        testTrace =
                                std::abs(trace - 1) < tol;
                        BOOST_CHECK(testTrace);

                        //check that this vector is normal to the edge \hat{e} adjacent by a to e
                        const int node = gel->SideNodeLocIndex(edgeIndex,iVec);
                        int neighEdgeIndex = -1;
                        for(int iNeigh = 0; iNeigh < nEdges; iNeigh++){
                            const int candidateEdgeIndex = iNeigh + nNodes;
                            for(int iNode = 0; iNode < 2; iNode++){
                                const bool cond1 = gel->SideNodeLocIndex(candidateEdgeIndex,iNode) == node;
                                const bool cond2 = candidateEdgeIndex != edgeIndex;
                                if(cond1 && cond2){
                                    neighEdgeIndex = candidateEdgeIndex;
                                }
                            }
                        }
                        if(neighEdgeIndex == -1){
                            DebugStop();
                        }
                        const int iNeighEdge = neighEdgeIndex - nNodes;

                        TPZVec<REAL> neighEdgeTgVector(3, 0);
                        for (auto i = 0; i < 3; i++) neighEdgeTgVector[i] = edgeTangentVecs(iNeighEdge,i);
                        const REAL neighEdgeLength = edgeLengthVec[iNeighEdge];
                        REAL neighTrace = 0;
                        for(auto x = 0; x < 3; x++) neighTrace += neighEdgeTgVector[x] * deformedDirections(x,2*iEdge+iVec);
                        neighTrace *=neighEdgeLength;
                        bool testNeighTrace = std::abs(neighTrace) < tol;
                        if(!testNeighTrace){
#ifdef NOISY_HCURL
                            std::cout<<"\t\t\tERROR"<<std::endl;
                            std::cout<<"\t\t\tintegral of tangential trace over neighbour edge: "<<neighTrace<<std::endl;
                            std::cout<<"\t\t\tneighbour tangent vector:"<<std::endl;
                            for (auto x = 0; x < 3; x++) std::cout << "\t\t" << neighEdgeTgVector[x]<< "\t";
                            std::cout<<std::endl;
#endif
                        }
                        BOOST_CHECK(testNeighTrace);
                    }//iVec

                    //testing neighbour elements
                    auto edgeType = TTopol::Type(edgeIndex);

                    const auto nEdgeNodes = MElementType_NNodes(edgeType);
                    TPZManVector<int64_t,4>edgeNodes(nEdgeNodes,-1);
                    for(auto iNode = 0; iNode < nEdgeNodes; iNode++) edgeNodes[iNode] = gel->SideNodeIndex(edgeIndex,iNode);

                    TPZGeoEl *edgeGel = gmesh->CreateGeoElement(edgeType, edgeNodes, 1, index, 0);
                    TPZFMatrix<REAL> edgeDeformedDirections(0,0),edgeMasterDirections(0,0);
                    TPZManVector<int,21> edgeTransformationIds(nSides - nNodes,pztopology::TPZLine::GetTransformId(nEdgeNodes, edgeNodes));
                    edgeTransformationIds[0] = pztopology::TPZLine::GetTransformId(nEdgeNodes, edgeNodes);
                    auxiliaryfuncs::ComputeDirections(edgeGel,edgeDeformedDirections,edgeMasterDirections,edgeTransformationIds);
                    for(auto iVec = 0; iVec <  2; iVec++) {// 2 v^{e,a}
                        REAL elTrace{0}, neighTrace{0};
                        for (auto x = 0; x < 3; x++) elTrace += edgeTgVector[x] * deformedDirections(x, 2 * iEdge + iVec);
                        for (auto x = 0; x < 3; x++) neighTrace += edgeTgVector[x] * edgeDeformedDirections(x, iVec);
                        const auto diffTrace = elTrace - neighTrace;
                        testTrace = std::abs(diffTrace) < tol;
                        BOOST_CHECK(testTrace);
                    }
                    for(auto iVec = 0; iVec <  1; iVec++) {// 1 v^{e,T}
                        REAL elTrace{0}, neighTrace{0};
                        for (auto x = 0; x < 3; x++) elTrace += edgeTgVector[x] * deformedDirections(x, 2 * nEdges + iEdge);
                        for (auto x = 0; x < 3; x++) neighTrace += edgeTgVector[x] * edgeDeformedDirections(x, 2);
                        const auto diffTrace = elTrace - neighTrace;
                        testTrace = std::abs(diffTrace) < tol;
                        BOOST_CHECK(testTrace);
                    }
                }//iEdge

                //calculating indexes and stuff
                TPZManVector<int> firstVfeVec(nFaces,-1);
                TPZManVector<TPZStack<int>> faceEdges(nFaces,TPZStack<int>(0,0));

                {
                    const int nEdgeVectors = nEdges * 3;
                    firstVfeVec[0] = nEdgeVectors;
                    for(auto iFace = 1; iFace < nFaces; iFace++){
                        TTopol::LowerDimensionSides(iFace - 1 + nEdges + nNodes, faceEdges[iFace-1], 1);
                        const int nFaceEdges = faceEdges[iFace-1].size();
                        firstVfeVec[iFace] = firstVfeVec[iFace - 1] + nFaceEdges;
                    }
                    TTopol::LowerDimensionSides(nFaces - 1 + nEdges + nNodes, faceEdges[nFaces-1], 1);
                }
                const int firstVftVec = firstVfeVec[nFaces-1] + faceEdges[nFaces-1].size();
                const int firstVfOrthVec = firstVftVec + 2 * nFaces;


                for(auto iFace = 0; iFace < nFaces; iFace++){
                    const int faceIndex = iFace + nEdges + nNodes;
                    auto faceType = TTopol::Type(faceIndex);
                    std::cout<<"\tface "<<iFace+1<<" out of "<<nFaces<<std::endl;
                    std::cout<<"\tface index "<<faceIndex<<std::endl;
                    std::cout<<"\tface type "<<MElementType_Name(faceType)<<std::endl;

                    const auto nFaceNodes = MElementType_NNodes(faceType);
                    TPZManVector<int64_t,4>faceNodes(nFaceNodes,-1);
                    for(auto iNode = 0; iNode < nFaceNodes; iNode++) faceNodes[iNode] = gel->SideNodeIndex(faceIndex,iNode);

                    TPZGeoEl *faceGel = gmesh->CreateGeoElement(faceType, faceNodes, 1, index, 0);
                    const auto nFaceEdges = faceGel->NSides() - nFaceNodes - 1;
                    TPZFMatrix<REAL> faceDeformedDirections(0,0),faceMasterDirections(0,0);
                    TPZManVector<int,21> faceTransformationIds(nSides - nNodes,-1);
                    for(auto iSide = 0 ; iSide < nFaceEdges+1; iSide++){
                        switch(faceType){
                            case ETriangle:
                                faceTransformationIds[iSide] = pztopology::TPZTriangle::GetTransformId(nFaceNodes + iSide, faceNodes);
                                break;
                            case EQuadrilateral:
                                faceTransformationIds[iSide] = pztopology::TPZQuadrilateral::GetTransformId(nFaceNodes + iSide, faceNodes);
                                break;
                            default:
                                DebugStop();
                        }
                    }
                    auxiliaryfuncs::ComputeDirections(faceGel,faceDeformedDirections,faceMasterDirections,faceTransformationIds);

                    //compute normal
                    TPZManVector<REAL,3> faceNormal(3,0);
                    auxiliaryfuncs::ComputeNormal(gmesh,faceNodes,faceNormal);


                    std::cout<<"\ttesting vfe vectors"<<std::endl;
                    const int nVfe = faceEdges[iFace].size();
                    const int firstVfeEl = firstVfeVec[iFace];
                    const int firstVfeFace = 3 * nVfe;
                    for(int iVec = 0; iVec < nVfe; iVec++){
#ifdef NOISY_HCURL
                        std::cout<<"\t\tvec "<<iVec+1<<" out of "<<nVfe<<std::endl;
#endif
                        TPZManVector<REAL,3> vfeElement(3,-1), vfeFace(3,-1);
                        for(auto x = 0; x < 3; x++) vfeElement[x] = deformedDirections(x,firstVfeEl + iVec);
                        for(auto x = 0; x < 3; x++) vfeFace[x] = faceDeformedDirections(x,firstVfeFace + iVec);


                        TPZManVector<REAL,3> temp(3,0);
                        TPZManVector<REAL,3> elTrace(3,0),faceTrace(3,0);
                        auxiliaryfuncs::VectorProduct(faceNormal,vfeElement,temp);
                        auxiliaryfuncs::VectorProduct(faceNormal,temp,elTrace);

                        auxiliaryfuncs::VectorProduct(faceNormal,vfeFace,temp);
                        auxiliaryfuncs::VectorProduct(faceNormal,temp,faceTrace);
                        REAL diff = 0;
                        for(auto x = 0; x < 3; x++) diff += (faceTrace[x]-elTrace[x])*(faceTrace[x]-elTrace[x]);
                        diff = sqrt(diff);
                        bool checkTraces = diff < tol;

                        if(!checkTraces){
                            std::cout<<"****************************************"<<std::endl;
                            std::cout<<"ERROR IN VFE VECTORS"<<std::endl;
                            std::cout<<"EL TYPE "<< MElementType_Name(elType)<<std::endl;
                            std::cout<<"FACE INDEX "<<faceIndex<<std::endl;
                            std::cout<<"FACE TYPE "<<MElementType_Name(faceType)<<std::endl;
                            std::cout<<"EDGE INDEX "<<faceEdges[iFace][iVec]<<std::endl;
                            std::cout<<"\t\tVECTOR (3D elem) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<vfeElement[x];
                            std::cout<<std::endl;
                            std::cout<<"\t\tVECTOR (2D face) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<vfeFace[x];
                            std::cout<<std::endl;
                            std::cout<<"\t\tTRACE (3D elem) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<elTrace[x];
                            std::cout<<std::endl;
                            std::cout<<"\t\tTRACE (2D face) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<faceTrace[x];
                            std::cout<<std::endl;
                            std::cout<<"****************************************"<<std::endl;
                        }
                        BOOST_CHECK(checkTraces);
                    }
                    std::cout<<"\ttesting vft vectors"<<std::endl;
                    const int nVft = 2;
                    const int firstVftEl = firstVftVec + 2 * iFace;
                    const int firstVftFace = 3 * nVfe + nVfe;
                    for(int iVec = 0; iVec < nVft; iVec++){
#ifdef NOISY_HCURL
                        std::cout<<"\t\tvec "<<iVec+1<<" out of "<<nVft<<std::endl;
#endif
                        TPZManVector<REAL,3> vftElement(3,-1), vftFace(3,-1);
                        for(auto x = 0; x < 3; x++) vftElement[x] = deformedDirections(x,firstVftEl + iVec);
                        for(auto x = 0; x < 3; x++) vftFace[x] = faceDeformedDirections(x,firstVftFace + iVec);


                        TPZManVector<REAL,3> temp(3,0);
                        TPZManVector<REAL,3> elTrace(3,0),faceTrace(3,0);
                        auxiliaryfuncs::VectorProduct(faceNormal,vftElement,temp);
                        auxiliaryfuncs::VectorProduct(faceNormal,temp,elTrace);

                        auxiliaryfuncs::VectorProduct(faceNormal,vftFace,temp);
                        auxiliaryfuncs::VectorProduct(faceNormal,temp,faceTrace);
                        REAL diff = 0;
                        for(auto x = 0; x < 3; x++) diff += (faceTrace[x]-elTrace[x])*(faceTrace[x]-elTrace[x]);
                        diff = sqrt(diff);
                        bool checkTraces = diff < tol;

                        if(!checkTraces){
                            std::cout<<"****************************************"<<std::endl;
                            std::cout<<"ERROR IN VFT VECTORS"<<std::endl;
                            std::cout<<"EL TYPE "<< MElementType_Name(elType)<<std::endl;
                            std::cout<<"FACE INDEX "<<faceIndex<<std::endl;
                            std::cout<<"FACE TYPE "<<MElementType_Name(faceType)<<std::endl;
                            std::cout<<"\t\tVECTOR (3D elem) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<vftElement[x];
                            std::cout<<std::endl;
                            std::cout<<"\t\tVECTOR (2D face) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<vftFace[x];
                            std::cout<<std::endl;
                            std::cout<<"\t\tTRACE (3D elem) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<elTrace[x];
                            std::cout<<std::endl;
                            std::cout<<"\t\tTRACE (2D face) =";
                            for(auto x = 0; x < 3; x++)std::cout<<"\t"<<faceTrace[x];
                            std::cout<<std::endl;
                            std::cout<<"****************************************"<<std::endl;
                        }
                        BOOST_CHECK(checkTraces);
                    }

                }


                delete gmesh;
            }//iPermute
        }//hcurltest::CompareVectorTraces

        void TestTracesUniformMesh(MElementType type, const int pOrder, const int dim){
            static std::string testName = __PRETTY_FUNCTION__;
            std::cout << testName << std::endl;
            std::cout<<"\t"<<MElementType_Name(type)<<" p = "<<pOrder<<std::endl;
            constexpr int ndiv{1};
            TPZManVector<int,2> matIds(2,-1);
            TPZGeoMesh *gmesh = [&]() -> TPZGeoMesh* {
                switch(dim){
                    case 2:
                        return auxiliaryfuncs::CreateGeoMesh2D(ndiv,ndiv,type,matIds);
                    case 3:
                        return auxiliaryfuncs::CreateGeoMesh3D(ndiv,ndiv,ndiv,type,matIds);
                    default:
                        DebugStop();
                }
                return nullptr;
            }();

#ifdef NOISY_HCURL
            std::ofstream outTXT("gmesh_"+MElementType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
            gmesh->Print(outTXT);
            outTXT.close();
#endif
#ifdef NOISY_HCURL_VTK
            std::ofstream outVTK("gmesh_"+MElementType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK);
            outVTK.close();
#endif
            auto cmesh = new TPZCompMesh(gmesh);
            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetDimModel(dim);

            auto mat = new TPZMatHelmholtz(dim,matIds[0],1,1);
            cmesh->InsertMaterialObject(mat);
            auto bcond = mat->CreateBC(mat,matIds[1],0,TPZFNMatrix<1,STATE>(1,1,0),TPZFNMatrix<1,STATE>(1,1,0));
            cmesh->InsertMaterialObject(bcond);
            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();

#ifdef NOISY_HCURL
            {
                std::ofstream outTXT("cmesh_"+MElementType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
                cmesh->Print(outTXT);
                outTXT.close();
            }
#endif

            ///this lambda expression will help calculating the tangential traces
            auto CheckTracesFunc = [](REAL &diffTrace, const TPZVec<REAL> &elShapeFunc, const TPZVec<REAL>&neighShapeFunc,
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
                        auxiliaryfuncs::VectorProduct(vec,elShapeFunc,temp);
                        auxiliaryfuncs::VectorProduct(vec,temp,elTrace);

                        auxiliaryfuncs::VectorProduct(vec,neighShapeFunc,temp);
                        auxiliaryfuncs::VectorProduct(vec,temp,neighTrace);
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
            };


            for(auto dummyCel : cmesh->ElementVec()){
                auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
                //skips boundary els
                if(!cel || cel->Reference()->Type() != type) continue;
                const int nState = cel->Material()->NStateVariables();
                const int elNNodes = MElementType_NNodes(type);
                for (auto iCon = 0; iCon <cel->NConnects(); iCon++) {
                    auto &con = cel->Connect(iCon);
                    const int nShape = con.NShape();
                    {
                        const bool check = con.NDof( *cmesh) == nShape * nState;
                        BOOST_CHECK_MESSAGE(check,"\n"+testName+" failed"+
                                                  "\ntopology: "+MElementType_Name(type)+"\n"
                        );
                    }

                    if(con.NElConnected() < 2) continue;

                    const int iSide = iCon + elNNodes;
                    const int pOrderIntRule = cel->EffectiveSideOrder(iSide)*2;
                    const auto gel = cel->Reference();
                    TPZIntPoints *sideIntRule = gel->CreateSideIntegrationRule(iSide, pOrderIntRule);
                    const int npts = sideIntRule->NPoints();

                    TPZTransform<> elTransform(gel->SideToSideTransform(iSide, gel->NSides() - 1));
                    TPZGeoElSide gelSide(gel, iSide);
                    const auto sideDim = gelSide.Dimension();


                    TPZGeoElSide neighGelSide = gelSide.Neighbour();
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
                    const int firstElShape = [&](){
                        int firstElShapeTemp = 0;
                        for(auto jCon = 0; jCon < iCon; jCon++){
                            firstElShapeTemp += cel->NConnectShapeF(jCon,cel->EffectiveSideOrder(jCon+elNNodes));
                        }
                        return firstElShapeTemp;
                    }();
                    const int nShapes = cel->NConnectShapeF(iCon,pOrder);
                    bool firstNeighbour{true};
//                    if(type == ETetraedro && gel->Id() == 0 && iSide == 5){
//                        std::cout<<"oi"<<std::endl;
//                    }
                    while(neighGelSide != gelSide) {
                        const auto neighCel = dynamic_cast<TPZInterpolatedElement *> (neighGelSide.Element()->Reference());
                        if (!neighCel) {
                            neighGelSide = neighGelSide.Neighbour();
                            continue;
                        }
                        const int neighNNodes = MElementType_NNodes(neighCel->Reference()->Type());
                        const auto neighSide = neighGelSide.Side();
                        TPZTransform<> neighTransform(neighCel->Reference()->SideToSideTransform(neighSide,
                                                                                                 neighCel->Reference()->NSides() -
                                                                                                 1));

                        {
                            const auto elConIndex = cel->SideConnectIndex(0, iSide);
                            const auto neighConIndex = neighCel->SideConnectIndex(0, neighSide);
                            const bool check = elConIndex == neighConIndex;
                            BOOST_CHECK_MESSAGE(check, "\n" + testName + " failed" +
                                                       "\ntopology: " + MElementType_Name(type) + "\n"
                            );
                        }
                        TPZMaterialData elData,neighData;
                        cel->InitMaterialData(elData);
                        neighCel->InitMaterialData(neighData);
                        TPZTransform<> localTransf(sideDim);
                        gelSide.SideTransform3(neighGelSide,localTransf);
                        const int firstNeighShape = [&](){
                            int firstNeighShapeTemp = 0;
                            const int neighCon = neighSide - neighGelSide.Element()->NNodes();
                            for(auto jCon = 0; jCon < neighCon; jCon++){
                                firstNeighShapeTemp += neighCel->NConnectShapeF(jCon,neighCel->EffectiveSideOrder(jCon+neighNNodes));
                            }
                            return firstNeighShapeTemp;
                        }();

                        const auto neighDim = neighGelSide.Element()->Dimension();
                        TPZManVector <REAL,3> pts(sideDim,0),ptsN(sideDim,0), ptEl(dim,0),ptNeigh(neighDim,0);
                        REAL w;
                        TPZFNMatrix<60,REAL> elShape,neighShape;
                        for (auto ipt = 0; ipt < npts; ipt++) {
                            sideIntRule->Point(ipt, pts, w);

                            elTransform.Apply(pts,ptEl);
                            cel->ComputeRequiredData(elData, ptEl);
                            TPZHCurlAuxClass::ComputeShape(elData.fVecShapeIndex, elData.phi,
                                                           elData.fDeformedDirections,elShape);
                            if(type == ETetraedro && pOrder ==4 && gel->Index()==4){
                                elData.jacobian.Print(std::cout);
                                elData.axes.Print(std::cout);
                                TPZFMatrix<REAL>curlPhi;
                                TPZHCurlAuxClass::ComputeCurl<3>(elData.fVecShapeIndex,elData.dphi,elData.fMasterDirections,elData.jacobian,elData.detjac,elData.axes,curlPhi);
                            }
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

                                BOOST_CHECK_MESSAGE(checkPhis,"\n"+testName+" failed: phis are different!"+
                                                                "\ntopology: "+MElementType_Name(type)+"\n"+
                                                                "side: "+std::to_string(iSide)+"\n"+
                                                                "p order: "+std::to_string(pOrder)+"\n"+
                                                                "elem phi index: "+std::to_string(elPhiIndex)+"\n"+
                                                                "neig phi index: "+std::to_string(neighPhiIndex)+"\n"
                                );
                                anyWrongCheck = !checkPhis || anyWrongCheck;
                                if(anyWrongCheck) {
                                    break;
                                }
                                for (auto x = 0; x < 3; x++) elShapeFunc[x] = elShape(elPhiIndex,x);
                                for (auto x = 0; x < 3; x++) neighShapeFunc[x] = neighShape(neighPhiIndex,x);
                                REAL diffTrace{0};
                                TPZManVector<REAL,3> elTrace,neighTrace;
                                const bool checkTraces = CheckTracesFunc(diffTrace,elShapeFunc,neighShapeFunc,vec,sideDim,elTrace,neighTrace);

                                {
                                    std::ostringstream traceMsg;
                                    if(!checkTraces){
                                        traceMsg <<"el func: ";
                                        for (auto x = 0; x < 3; x++) traceMsg << elShapeFunc[x]<<"\t";
                                        traceMsg<<"\nneigh func: ";
                                        for (auto x = 0; x < 3; x++) traceMsg << neighShapeFunc[x]<<"\t";
                                        traceMsg <<"\n";
                                        traceMsg <<"el trace: ";
                                        for (auto x = 0; x < elTrace.size(); x++) traceMsg << elTrace[x]<<"\t";
                                        traceMsg<<"\nneigh trace: ";
                                        for (auto x = 0; x < neighTrace.size(); x++) traceMsg << neighTrace[x]<<"\t";
                                        traceMsg <<"\n";
                                    }
                                    BOOST_CHECK_MESSAGE(checkTraces,"\n"+testName+" failed"+
                                                                    "\ntopology: "+MElementType_Name(type)+"\n"+
                                                                    "side: "+std::to_string(iSide)+"\n"+
                                                                    "p order: "+std::to_string(pOrder)+"\n"+
                                                                    "diff traces: "+std::to_string(diffTrace)+"\n"+
                                                                    traceMsg.str()
                                    );
                                }
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
                                    TPZStack<int> subSides;
                                    gel->LowerDimensionSides(gelSide.Side(),subSides);
                                    for(auto iSubSide = gelSide.NSideNodes(); iSubSide < gelSide.NSides() - 1; iSubSide++){
                                        const auto subSide = subSides[iSubSide];
                                        firstSideShapeTemp += cel->NConnectShapeF(subSide - gel->NNodes(),pOrder);
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
                                    {
                                        std::ostringstream traceMsg;
                                        if(!checkTraces){
                                            traceMsg <<"el func: ";
                                            for (auto x = 0; x < 3; x++) traceMsg << elShapeFunc[x]<<"\t";
                                            traceMsg<<"\nside func: ";
                                            for (auto x = 0; x < 3; x++) traceMsg << sideShapeFunc[x]<<"\t";
                                            traceMsg <<"\n";
                                            traceMsg <<"el trace: ";
                                            for (auto x = 0; x < elTrace.size(); x++) traceMsg << elTrace[x]<<"\t";
                                            traceMsg<<"\nside trace: ";
                                            for (auto x = 0; x < elTrace.size(); x++) traceMsg << sideTrace[x]<<"\t";
                                            traceMsg <<"\n";
                                        }
                                        BOOST_CHECK_MESSAGE(checkTraces,"\n"+testName+" failed"+
                                                                        "\ntopology: "+MElementType_Name(type)+"\n"+
                                                                        "side: "+std::to_string(iSide)+"\n"+
                                                                        "p order: "+std::to_string(pOrder)+"\n"+
                                                                        "diff traces: "+std::to_string(diffTrace)+"\n"+
                                                                        traceMsg.str()
                                        );
                                    }
                                    anyWrongCheck = !checkTraces || anyWrongCheck;
                                }
                                if(anyWrongCheck){
                                    break;
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
//            TPZAnalysis an(cmesh,false);
//            const int postProcessResolution = 3;
//            const std::string executionInfo = [&](){
//                std::string name("");
//                name.append(MElementType_Name(type));
//                name.append(std::to_string(pOrder));
//                return name;
//            }();
//
//            const std::string plotfile = "solution"+executionInfo+".vtk";//where to print the vtk files
//            TPZStack<std::string> scalnames, vecnames;
//            vecnames.Push("E");//print the state variable
//            auto sol = an.Solution();
//            sol.Zero();
//            for(int i = 0; i < sol.Rows(); i++){
//                sol(i - 1 < 0 ? 0 : i - 1 , 0) = 0;
//                sol(i,0) = 1;
//                an.LoadSolution(sol);
//                an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
//                an.SetStep(i);
//                an.PostProcess(postProcessResolution);
//            }
            delete cmesh;
            delete gmesh;
        }//hcurltest::TestTracesUniformMesh


        template <class TGEOM>
        void PrintShapeFunctions(const int pOrder){
            constexpr auto nNodes{TGEOM::NCornerNodes};
            constexpr auto dim{TGEOM::Dimension};
            TPZGeoMesh *gmesh = new TPZGeoMesh();
            gmesh->SetDimension(dim);
            //Auxiliary vector to store coordinates:
            TPZManVector<REAL,3> coords(3, 0.);
            for(int iNode = 0 ; iNode < nNodes; iNode++){
                TGEOM::ParametricDomainNodeCoord(iNode,coords);
                auto newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coords, *gmesh);
            }
            TPZManVector<int64_t,TGEOM::NCornerNodes> nodesIdVec(TGEOM::NCornerNodes,-1);
            for(int i = 0 ; i < TGEOM::NCornerNodes; i++) nodesIdVec[i] = i;
            int64_t index{-1};
            TPZGeoEl *gel = gmesh->CreateGeoElement(TGEOM::Type(), nodesIdVec, 1, index, 0);
            gmesh->BuildConnectivity();
            TPZCompMesh *cmesh = new TPZCompMesh (gmesh);

            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetDimModel(dim);

            auto mat = new TPZMatHelmholtz(dim,1,1,1);
            cmesh->InsertMaterialObject(mat);
            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->CleanUpUnconnectedNodes();

            TPZAnalysis an(cmesh,false);
            const int postProcessResolution = 3;
            const std::string executionInfo = [&](){
                std::string name("");
                name.append(MElementType_Name(TGEOM::Type()));
                name.append(std::to_string(pOrder));
                return name;
            }();

            const std::string plotfile = "shapeFuncs"+executionInfo+".vtk";//where to print the vtk files
            TPZStack<std::string> scalnames, vecnames;
            vecnames.Push("E");//print the state variable
            auto sol = an.Solution();
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


        namespace auxiliaryfuncs{
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

#ifdef NOISY_HCURL
                std::cout << std::endl;
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

            void ComputeNormal(TPZGeoMesh *gmesh, TPZVec<int64_t> faceNodes,TPZVec<REAL> &result){
                TPZVec<REAL> v1(3);
                TPZVec<REAL> v2(3);
                TPZVec<REAL> normal(3);
                TPZVec<REAL> p0(3),p1(3),p2(3);
                gmesh->NodeVec()[faceNodes[0]].GetCoordinates(p0);
                gmesh->NodeVec()[faceNodes[1]].GetCoordinates(p1);
                gmesh->NodeVec()[faceNodes[2]].GetCoordinates(p2);
                v1[0]=p1[0]-p0[0];
                v1[1]=p1[1]-p0[1];
                v1[2]=p1[2]-p0[2];
                v2[0]=p2[0]-p1[0];
                v2[1]=p2[1]-p1[1];
                v2[2]=p2[2]-p1[2];
                VectorProduct(v1,v2,result);
                REAL norm = 0;
                for(auto x = 0; x < 3; x++) norm += result[x] * result[x];
                norm = std::sqrt(norm);
                for(auto x = 0; x < 3; x++) result[x] /= norm;
            };//ComputeNormal

            TPZGeoMesh *CreateGeoMesh2D(int nelx, int nely, MElementType meshType, TPZVec<int> &matIds) {
                //Creating geometric mesh, nodes and elements.
                //Including nodes and elements in the mesh object:
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                constexpr int dim{2};
                gmesh->SetDimension(dim);

                //Auxiliary vector to store coordinates:
                TPZVec<REAL> coord1(3, 0.);
                TPZVec<REAL> coord2(3, 0.);
                coord1[0] = 0;coord1[1] = 0;coord1[2] = 0;
                coord2[0] = 1;coord2[1] = 1;coord2[2] = 0;

                TPZManVector<int> nelem(2, 1);
                nelem[0] = nelx;
                nelem[1] = nely;

                TPZGenGrid gengrid(nelem, coord1, coord2);

                switch (meshType) {
                    case EQuadrilateral:
                        gengrid.SetElementType(EQuadrilateral);
                        break;
                    case ETriangle:
                        gengrid.SetElementType(ETriangle);
                        break;
                    default:
                        DebugStop();
                }
                constexpr int matIdDomain = 1, matIdBoundary = 2;
                gengrid.Read(gmesh, matIdDomain);
                gengrid.SetBC(gmesh, 4, matIdBoundary);
                gengrid.SetBC(gmesh, 5, matIdBoundary);
                gengrid.SetBC(gmesh, 6, matIdBoundary);
                gengrid.SetBC(gmesh, 7, matIdBoundary);

                gmesh->BuildConnectivity();

                {
                    TPZCheckGeom check(gmesh);
                    check.CheckUniqueId();
                }
                matIds.Resize(2);
                matIds[0] = matIdDomain;
                matIds[1] = matIdBoundary;
                //Printing geometric mesh:

                //ofstream bf("before.vtk");
                //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
                return gmesh;
            }//CreateGeoMesh2D

            TPZGeoMesh *CreateGeoMesh3D(int nelx, int nely, int nelz, MElementType meshType, TPZVec<int> &matIds) {
                //Creating geometric mesh, nodes and elements.
                //Including nodes and elements in the mesh object:
                //create boundary elements
                constexpr int minX{0},minY{0},minZ{0};
                constexpr int maxX{1},maxY{1},maxZ{1};
                constexpr int matIdDomain{1};
                constexpr int matIdBoundary{2};

                TPZGenGrid3D genGrid3D(minX,minY,minZ,maxX,maxY,maxZ,nelx,nely,nelz,meshType);
                genGrid3D.BuildVolumetricElements(matIdDomain);
                TPZGeoMesh *gmesh = genGrid3D.BuildBoundaryElements(matIdBoundary,matIdBoundary,matIdBoundary,matIdBoundary,matIdBoundary,matIdBoundary);
                gmesh->BuildConnectivity();

                {
                    TPZCheckGeom check(gmesh);
                    check.CheckUniqueId();
                    check.PerformCheck();
                }
                matIds.Resize(2);
                matIds[0] = matIdDomain;
                matIds[1] = matIdBoundary;
                //Printing geometric mesh:

                //ofstream bf("before.vtk");
                //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
                return gmesh;
            }//CreateGeoMesh3D
        }//namespace
    }//namespace


BOOST_AUTO_TEST_SUITE_END()

#ifdef NOISY_HCURL
#undef NOISY_HCURL
#endif
#ifdef NOISYVTK_HCURL
#undef NOISYVTK_HCURL
#endif

#endif

