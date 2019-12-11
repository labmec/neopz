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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testhcurl"));
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

#include <pzgengrid.h>
#include <pzcmesh.h>
#include <TPZMatHelmholtz2D.cpp>
#include <pzanalysis.h>
#include <pzintel.h>
#include <TPZCompElHCurl.h>
#include <pzgeotetrahedra.h>
#include <pzgeoprism.h>
#include <TPZGeoCube.h>
#include <pzgeoelrefless.h>
#include <TPZVTKGeoMesh.h>
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
#endif

#define NOISY_HCURL //outputs useful debug info
#define NOISY_HCURL_VTK
//std::string dirname = PZSOURCEDIR;

#ifdef USING_BOOST
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

        namespace auxiliaryfuncs{
            void ComputeDirections (TPZGeoEl *, TPZFMatrix<REAL> &, TPZFMatrix<REAL> &,const TPZVec<int> &);

            void VectorProduct(TPZVec<REAL> &, TPZVec<REAL> &, TPZVec<REAL> &);

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

            auto mat = new TPZMatHelmholtz2D(matIds[0],1,1);
            cmesh->InsertMaterialObject(mat);

            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();
            for(auto dummyCel : cmesh->ElementVec()){
                auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
                if(!cel) continue;
                int nState = cel->Material()->NStateVariables();
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

                    const int iSide = iCon + MElementType_NNodes(type);
                    const int pOrderIntRule = [&](){
                        auto dummy = dynamic_cast<TPZInterpolatedElement *>(cel);
                        auto pOrder = dummy->EffectiveSideOrder(iSide)*2;
                        return pOrder;
                    }();
                    TPZIntPoints *sideIntRule = cel->Reference()->CreateSideIntegrationRule(iSide, pOrderIntRule);
                    TPZTransform<> elTransform(cel->Reference()->SideToSideTransform(iSide, cel->Reference()->NSides() - 1));
                    const int npts = sideIntRule->NPoints();
                    TPZGeoElSide gelSide(cel->Reference(), iSide);
                    TPZGeoElSide neighGelSide = gelSide.Neighbour();
                    auto sideDim = gelSide.Dimension();
                    auto gel = gelSide.Element();
                    //the following vector will be the edge tg vector if 2D, the normal vector if 3D
                    TPZManVector<REAL,3> vec(3,0);
                    switch(sideDim){
                        case 0: continue;
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
                    int firstElShape = 0;
                    for(auto jCon = 0; jCon < iCon; jCon++)   firstElShape += cel->NConnectShapeF(jCon,pOrder);
                    const int nShapes = cel->NConnectShapeF(iCon,pOrder);

                    while(neighGelSide != gelSide && sideDim != 0){
                        auto neighCel = dynamic_cast<TPZInterpolatedElement *> (neighGelSide.Element()->Reference());
                        auto neighSide = neighGelSide.Side();
                        TPZTransform<> neighTransform(neighCel->Reference()->SideToSideTransform(neighSide,
                                neighCel->Reference()->NSides() - 1));
                        TPZMaterialData elData,neighData;
                        cel->InitMaterialData(elData);
                        neighCel->InitMaterialData(neighData);
                        TPZTransform<> localTransf(sideDim);
                        neighGelSide.SideTransform3(gelSide,localTransf);
                        TPZManVector <REAL,3> pts(sideDim,0),ptsN(sideDim,0), ptEl(dim,0),ptNeigh(dim,0);
                        TPZFNMatrix<30,STATE> elShape,neighShape;

                        int firstNeighShape = 0;
                        const int neighCon = neighSide - neighGelSide.Element()->NNodes();
                        for(auto jCon = 0; jCon < neighCon; jCon++)   firstNeighShape += cel->NConnectShapeF(jCon,pOrder);

                        for (auto ipt = 0; ipt < npts; ipt++) {
                            REAL w;

                            sideIntRule->Point(ipt, pts, w);


                            elTransform.Apply(pts,ptEl);
                            cel->ComputeRequiredData(elData, ptEl);
                            cel->ComputeShape(ptEl,elData);
                            TPZHCurlAuxClass::ComputeShape(elData.fVecShapeIndex, elData.phi,
                                                           elData.fDeformedDirections,elShape);
                            localTransf.Apply(pts,ptsN);
                            neighTransform.Apply(ptsN,ptNeigh);
                            neighCel->ComputeRequiredData(neighData, ptNeigh);
                            neighCel->ComputeShape(ptNeigh,neighData);
                            TPZHCurlAuxClass::ComputeShape(neighData.fVecShapeIndex, neighData.phi,
                                                           neighData.fDeformedDirections,neighShape);

                            TPZManVector<REAL,3> elShapeFunc(3,0), neighShapeFunc(3,0);
                            for(auto iShape = 0; iShape < nShapes; iShape ++){
                                for (auto x = 0; x < 3; x++) elShapeFunc[x] = elShape(firstElShape + iShape,x);
                                for (auto x = 0; x < 3; x++) neighShapeFunc[x] = neighShape(firstNeighShape + iShape,x);
                                REAL elTrace{0}, neighTrace{0};
                                const bool checkTraces = [&](REAL &elTrace, REAL &neighTrace){
                                    switch(sideDim){
                                        case 1:{
                                            for (auto x = 0; x < 3; x++) elTrace += elShapeFunc[x]*vec[x];
                                            for (auto x = 0; x < 3; x++) neighTrace += neighShapeFunc[x]*vec[x];
                                            return std::abs(elTrace - neighTrace) < tol;
                                        }
                                            break;
                                        case 2:{
                                            TPZManVector<REAL,3> temp(3,0);
                                            TPZManVector<REAL,3> elTrace(3,0),faceTrace(3,0);

                                            auxiliaryfuncs::VectorProduct(vec,elShapeFunc,temp);
                                            auxiliaryfuncs::VectorProduct(vec,temp,elTrace);

                                            auxiliaryfuncs::VectorProduct(vec,neighShapeFunc,temp);
                                            auxiliaryfuncs::VectorProduct(vec,temp,faceTrace);
                                            REAL diff = 0;
                                            for(auto x = 0; x < 3; x++) diff += (faceTrace[x]-elTrace[x])*(faceTrace[x]-elTrace[x]);
                                            diff = sqrt(diff);
                                            return diff < tol;
                                        }
                                            break;
                                        default:
                                            return false;
                                    }
                                }(elTrace,neighTrace);
                                BOOST_CHECK_MESSAGE(checkTraces,"\n"+testName+" failed"+
                                                          "\ntopology: "+MElementType_Name(type)+"\n"+
                                                          "side: "+std::to_string(iSide)+"\n"+
                                                          "p order: "+std::to_string(pOrder)+"\n"
                                                          "el trace: "+std::to_string(elTrace)+"\n"
                                                          "neigh trace: "+std::to_string(neighTrace)+"\n"
                                );
                            }
                        }
                        neighGelSide = neighGelSide.Neighbour();
                    }
                }
            }
            /***********************************************************************************************************
             *              the following lines might be useful for analysing the basis functions
            ***********************************************************************************************************/
//            TPZAnalysis an(cmesh);
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

            void VectorProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result){
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
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                constexpr int dim{3};
                gmesh->SetDimension(dim);

                //Auxiliary vector to store coordinates:
                constexpr REAL maxX{1},maxY{1},maxZ{1};

                //create nodes
                [&](){
                    TPZManVector<REAL,3> nodeCoords(3,0);
                    int64_t newindex = -1;
                    for(auto iZ = 0; iZ <= nelz; iZ++){
                        for(auto iY = 0; iY <= nely; iY++){
                            for(auto iX = 0; iX <= nelx; iX++){
                                nodeCoords[0] = (REAL)iX/nelx;
                                nodeCoords[1] = (REAL)iY/nely;
                                nodeCoords[2] = (REAL)iZ/nelz;
                                newindex = gmesh->NodeVec().AllocateNewElement();
                                gmesh->NodeVec()[newindex].Initialize(nodeCoords, *gmesh);
                            }
                        }
                    }
                }();


                constexpr int matIdDomain = 1, matIdBoundary = 2;

                //create elements
                [&](){
                    TPZGeoEl * gel = nullptr;
                    TPZManVector<int64_t,8> nodesIdVec(MElementType_NNodes(meshType),-1);
                    for(auto iZ = 0; iZ < nelz; iZ++){
                        for(auto iY = 0; iY < nely; iY++){
                            for(auto iX = 0; iX < nelx; iX++){
                                const auto firstNodeId = iZ * (nelx+1) * (nely+1) + iY * (nelx+1) + iX;/*lower left node*/
                                switch (meshType) {
                                    case ETetraedro:
                                        nodesIdVec[0] = firstNodeId;
                                        nodesIdVec[1] = firstNodeId + 1;
                                        nodesIdVec[2] = firstNodeId + (nelx+1);
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1);
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*gmesh);
                                        nodesIdVec[0] = firstNodeId + 1;
                                        nodesIdVec[1] = firstNodeId + (nelx+1) * (nely+1) + 1;
                                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1);
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*gmesh);
                                        nodesIdVec[0] = firstNodeId + 1;
                                        nodesIdVec[1] = firstNodeId + (nelx+1) + 1;
                                        nodesIdVec[2] = firstNodeId + (nelx+1);
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*gmesh);
                                        nodesIdVec[0] = firstNodeId + (nelx+1);
                                        nodesIdVec[1] = firstNodeId + (nelx+1) * (nely+1);
                                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1);
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*gmesh);
                                        nodesIdVec[0] = firstNodeId + 1;
                                        nodesIdVec[1] = firstNodeId + (nelx+1);
                                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1);
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*gmesh);
                                        break;
                                    case ECube:
                                        nodesIdVec[0] = firstNodeId;
                                        nodesIdVec[1] = firstNodeId + 1;
                                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                                        nodesIdVec[3] = firstNodeId + (nelx+1);
                                        nodesIdVec[4] = firstNodeId + (nelx+1) * (nely+1);
                                        nodesIdVec[5] = firstNodeId + (nelx+1) * (nely+1) + 1;
                                        nodesIdVec[6] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        nodesIdVec[7] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1);
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoCube>(nodesIdVec,matIdDomain,*gmesh);
                                        break;
                                    case EPrisma:
                                        nodesIdVec[0] = firstNodeId;
                                        nodesIdVec[1] = firstNodeId + 1;
                                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1);
                                        nodesIdVec[4] = firstNodeId + (nelx+1) * (nely+1) + 1;
                                        nodesIdVec[5] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoPrism>(nodesIdVec,matIdDomain,*gmesh);

                                        nodesIdVec[0] = firstNodeId;
                                        nodesIdVec[1] = firstNodeId + (nelx+1);
                                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1);
                                        nodesIdVec[4] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1);
                                        nodesIdVec[5] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1) + 1;
                                        gel = new TPZGeoElRefLess<pzgeom::TPZGeoPrism>(nodesIdVec,matIdDomain,*gmesh);
                                        break;
                                    default:
                                        DebugStop();
                                }
                            }
                        }
                    }

                }();

                gmesh->BuildConnectivity();

                {
                    TPZCheckGeom check(gmesh);
                    check.CheckUniqueId();
                }
                matIds.Resize(1);
                matIds[0] = matIdDomain;
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

