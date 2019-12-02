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
// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz #define BOOST_TEST_MAIN pz hcurl_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"

#endif

//#define NOISY_HCURL //outputs useful debug info

//std::string dirname = PZSOURCEDIR;

#ifdef USING_BOOST

struct SuiteInitializer{
    SuiteInitializer(){
        InitializePZLOG();
    }
};
BOOST_FIXTURE_TEST_SUITE(hcurl_tests,SuiteInitializer)
    namespace hcurltest{
        constexpr REAL tol = 1e-10;
        template <class TTopol>
        void CompareVectorTraces(const TPZFMatrix<REAL> &);
        void TestExampleMesh2D(MElementType type);

        namespace auxiliaryfuncs{
            void ComputeDirections (TPZGeoEl *, TPZFMatrix<REAL> &, TPZFMatrix<REAL> &,const TPZVec<int> &);

            void VectorProduct(TPZVec<REAL> &, TPZVec<REAL> &, TPZVec<REAL> &);

            void ComputeNormal(TPZGeoMesh *, TPZVec<int64_t>, TPZVec<REAL> &);

            TPZGeoMesh *CreateGeoMesh(const int dim, int nelx, int nely, MElementType meshType, TPZVec<int> &matIds);
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
        hcurltest::TestExampleMesh2D(ETriangle);

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

        void TestExampleMesh2D(MElementType type){
            std::cout << __PRETTY_FUNCTION__ << MElementType_Name(type)<<std::endl;
            constexpr int dim{2};
            constexpr int ndiv{4};
            constexpr int pOrder{1};
            TPZManVector<int,2> matIds(2,-1);
            auto gmesh = auxiliaryfuncs::CreateGeoMesh(dim,ndiv,ndiv,type,matIds);
            auto cmesh = new TPZCompMesh(gmesh);
            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetDimModel(dim);

            auto mat = new TPZMatHelmholtz2D(matIds[0],1,1);
            cmesh->InsertMaterialObject(mat);
            //Boundary conditions
            constexpr int dirichlet = 0;
            constexpr int neumann = 1;
            TPZFMatrix<STATE> val1(1,1,0.0);
            TPZFMatrix<STATE> val2(1,1,0.0);

            const int &matIdBc1 = matIds[1];
            val2(0,0)=0.0;
            auto bc1 = mat->CreateBC(mat, matIdBc1, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc1);
            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();
            delete cmesh;
            delete gmesh;
        }//hcurltest::TestExampleMesh2D


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

            TPZGeoMesh *CreateGeoMesh(const int dim, int nelx, int nely, MElementType meshType, TPZVec<int> &matIds) {
                //Creating geometric mesh, nodes and elements.
                //Including nodes and elements in the mesh object:
                TPZGeoMesh *gmesh = new TPZGeoMesh();
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
            }//CreateGeoMesh
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

