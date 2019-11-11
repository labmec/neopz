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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testhcurl"));
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz blend_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"

#endif

#define NOISY_HCURL //outputs useful debug info
//#define NOISYVTK_HCURL//outputs even more debug info, printing relevant info in .vtk format

std::string dirname = PZSOURCEDIR;

#ifdef USING_BOOST


BOOST_AUTO_TEST_SUITE(hcurl_tests)
    namespace hcurltest{
        constexpr REAL tol = 1e-10;
        template <class TTopol>
        void ComparePermutedVectors(const TPZFMatrix<REAL> &);
        template <class TTopol>
        void CompareTraces(const TPZFMatrix<REAL> &);
    }

    BOOST_AUTO_TEST_CASE(hcurl_permutation_tests) {
        InitializePZLOG();
        {
            TPZFMatrix<REAL> nodeCoords(3,3);
            nodeCoords(0,0) = -1;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  0;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            hcurltest::ComparePermutedVectors<pztopology::TPZTriangle>(nodeCoords);
            hcurltest::CompareTraces<pztopology::TPZTriangle>(nodeCoords);
        }
        {
            TPZFMatrix<REAL> nodeCoords(4,3);
            nodeCoords(0,0) =  0;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  1;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            nodeCoords(3,0) =  0;   nodeCoords(3,1) =  1;   nodeCoords(3,2) =  0;
            hcurltest::ComparePermutedVectors<pztopology::TPZQuadrilateral>(nodeCoords);
            hcurltest::CompareTraces<pztopology::TPZQuadrilateral>(nodeCoords);
        }

        {
            TPZFMatrix<REAL> nodeCoords(4,3);
            nodeCoords(0,0) = -1;   nodeCoords(0,1) =  0;   nodeCoords(0,2) =  0;
            nodeCoords(1,0) =  1;   nodeCoords(1,1) =  0;   nodeCoords(1,2) =  0;
            nodeCoords(2,0) =  0;   nodeCoords(2,1) =  1;   nodeCoords(2,2) =  0;
            nodeCoords(3,0) =  0;   nodeCoords(3,1) =  0;   nodeCoords(3,2) =  1;
            hcurltest::ComparePermutedVectors<pztopology::TPZTetrahedron>(nodeCoords);
            hcurltest::CompareTraces<pztopology::TPZTetrahedron>(nodeCoords);
        }
    }

    namespace hcurltest{


        template <class TTopol>
        void ComparePermutedVectors(const TPZFMatrix<REAL> &nodeCoords) {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            auto elType = TTopol::Type();

            const int64_t nSides = TTopol::NumSides();
            const int64_t nNodes = TTopol::NumSides(0);
            const int64_t nEdges = TTopol::NumSides(1);
            const int64_t nFaces = TTopol::NumSides(2);
            const int64_t dim = TTopol::Dimension;

            //mesh used for calculating the transformation between the permuted element and the original one
            TPZGeoMesh *permuteGMesh = new TPZGeoMesh();
            //actual mesh
            TPZGeoMesh *gmesh = new TPZGeoMesh();
            ///creating the mesh nodes.
            TPZVec<REAL> nodeLocalCoord(0), xiNode(dim, 0), xNode(3, 0);
            int64_t newIndex = -1;
            for (int iNode = 0; iNode < nNodes; iNode++) {
                auto transf = TTopol::TransformSideToElement(iNode);
                transf.Apply(nodeLocalCoord, xiNode);
                newIndex = gmesh->NodeVec().AllocateNewElement();
                permuteGMesh->NodeVec().AllocateNewElement();
//#ifdef NOISY_HCURL
//                std::cout<<"node "<<iNode<<" coords:"<<std::endl;
//                for(auto x = 0; x < dim; x++) std::cout<<xiNode[x]<<"\t";
//                std::cout<<std::endl;
//#endif
                permuteGMesh->NodeVec()[newIndex].Initialize(xiNode, *permuteGMesh);
                for (auto x = 0; x < 3; x++) xNode[x] = nodeCoords.GetVal(iNode, x);
                gmesh->NodeVec()[newIndex].Initialize(xNode, *gmesh);
            }
            //create original element in dummy mesh
            TPZGeoEl *originalDummyEl = nullptr;
            {
                TPZVec<int64_t> nodesPerm(nNodes, 0);
                for (auto x = 0; x < nNodes; x++) nodesPerm[x] = x;
                int64_t index, matId = 1;
                originalDummyEl = permuteGMesh->CreateGeoElement(elType, nodesPerm, matId, index, 0);
            }

            const int64_t nPermutations = TTopol::NPermutations;

            auto CreatePermutedEl = [gmesh, permuteGMesh, elType, nSides, nNodes, originalDummyEl](
                    TPZVec<int> &currentPermutation, TPZFMatrix<REAL> &transMult) {
                TPZVec<int64_t> nodesPerm(nNodes, 0);
                for (auto x = 0; x < nNodes; x++) nodesPerm[x] = currentPermutation[x];
                int64_t index, matId = 1;
                TPZGeoEl *dummy = permuteGMesh->CreateGeoElement(elType, nodesPerm, matId, index, 0);
                TPZTransform<> transf = dummy->ComputeParamTrans(originalDummyEl, nSides - 1, nSides - 1);
//                TPZTransform<> transf = originalDummyEl->ComputeParamTrans(dummy,nSides-1,nSides-1);
                transMult = transf.Mult();
                return gmesh->CreateGeoElement(elType, nodesPerm, matId, index, 0);
            };

            int permute = 0;
            TPZVec<int> currentPermutation(nSides, 0);
            pztopology::GetPermutation<TTopol>(permute, currentPermutation);
            TPZFMatrix<REAL> transMult(3, dim, 0);
            TPZGeoEl *originalEl = CreatePermutedEl(currentPermutation, transMult);
            TPZFMatrix<REAL> gradX(3, dim, 0);
            for (auto x = 0; x < dim; x++) gradX(x, x) = 1;
            const int nVec = dim * nSides;
            TPZFMatrix<REAL> masterDirections(3, nVec), originalDirections(3, nVec);
            TTopol::ComputeHCurlDirections(gradX, masterDirections);//these are the vectors on the master element



            TPZVec<REAL> xiCenter(dim, 0);
            TPZFMatrix<REAL> originalJacobian(dim, dim, 0), originalAxes(dim, 3), originalJacInv(dim, dim, 0);
            REAL originalDetJac = 0;
            originalEl->JacobianXYZ(xiCenter, originalJacobian, originalAxes, originalDetJac, originalJacInv);

#ifdef NOISY_HCURL
            std::cout << std::endl;
            originalJacobian.Print("Original Jacobian:");
            originalAxes.Print("Original Axes:");
#endif
            for (auto iVec = 0; iVec < nVec; iVec++) {

                TPZManVector<REAL, 3> tempDirection(dim, 0);
                for (auto i = 0; i < dim; i++) {
                    //covariant piola transform: J^{-T}
                    tempDirection[i] = 0;
                    for (auto j = 0; j < dim; j++) tempDirection[i] += originalJacInv(j, i) * masterDirections(j, iVec);
                }
                for (auto i = 0; i < 3; i++) {
                    originalDirections(i, iVec) = 0;
                    for (auto j = 0; j < dim; j++) originalDirections(i, iVec) += originalAxes(j, i) * tempDirection[j];
                }
            }


            TPZGeoEl *permutedEl = nullptr;
            TPZFMatrix<REAL> permutedJacobian(dim, dim, 0), permutedAxes(dim, 3), permutedJacInv(dim, dim, 0);
            REAL permutedDetJac = 0;

            TTopol::CenterPoint(nSides - 1, xiCenter);
            TPZFMatrix<REAL> permutedDirections(3, nVec, 0);
            for (permute = 1; permute < nPermutations; permute++) {
                pztopology::GetPermutation<TTopol>(permute, currentPermutation);
                std::cout << "\tpermutation " << permute << " out of " << nPermutations << ". side ordering:"
                          << std::endl;
                for (auto iSide = 0; iSide < nSides; iSide++) std::cout << "\t" << currentPermutation[iSide];
                std::cout << std::endl;
                permutedEl = CreatePermutedEl(currentPermutation, transMult);

                //in all valid permutations the gradX is constant therefore the xi point should not matter
                permutedEl->JacobianXYZ(xiCenter, permutedJacobian, permutedAxes, permutedDetJac, permutedJacInv);

#ifdef NOISY_HCURL
                std::cout << std::endl;
                permutedJacobian.Print("Jacobian:");
                permutedAxes.Print("Axes:");
                transMult.Print("Permute Transform:");
#endif
                for (auto iVec = 0; iVec < nVec; iVec++) {

                    TPZManVector<REAL, 3> tempDirection(dim, 0), tempDirection2(dim, 0);
                    for (auto i = 0; i < dim; i++) {
                        //transpose of permutation matrix
                        tempDirection[i] = 0;
                        for (auto j = 0; j < dim; j++) tempDirection[i] += transMult(j, i) * masterDirections(j, iVec);
                    }
                    for (auto i = 0; i < dim; i++) {
                        //covariant piola transform: J^{-T}
                        tempDirection2[i] = 0;
                        for (auto j = 0; j < dim; j++) tempDirection2[i] += permutedJacInv(j, i) * tempDirection[j];
                    }
                    for (auto i = 0; i < 3; i++) {
                        permutedDirections(i, iVec) = 0;
                        for (auto j = 0; j < dim; j++)
                            permutedDirections(i, iVec) += permutedAxes(j, i) * tempDirection2[j];
                    }
                    REAL diff = 0;
                    for (auto x = 0; x < 3; x++)
                        diff += (originalDirections(x, iVec) - permutedDirections(x, iVec)) *
                                (originalDirections(x, iVec) - permutedDirections(x, iVec));
                    bool test = diff < tol;
                    BOOST_CHECK(test);
                }
            }
            delete permuteGMesh;
            delete gmesh;
        }

        template<class TTopol>
        void CompareTraces(const TPZFMatrix<REAL> &nodeCoords) {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            auto elType = TTopol::Type();

            const int64_t nSides = TTopol::NumSides();
            const int64_t nNodes = TTopol::NumSides(0);
            const int64_t nEdges = TTopol::NumSides(1);
            const int64_t nFaces = TTopol::NumSides(2);
            const int64_t dim = TTopol::Dimension;
            //actual mesh
            TPZGeoMesh *gmesh = new TPZGeoMesh();
            ///creating the mesh nodes.
            TPZVec<REAL> nodeLocalCoord(0), xNode(3, 0);
            int64_t newIndex = -1;
            for (int iNode = 0; iNode < nNodes; iNode++) {
                newIndex = gmesh->NodeVec().AllocateNewElement();
//#ifdef NOISY_HCURL
//                std::cout<<"node "<<iNode<<" coords:"<<std::endl;
//                for(auto x = 0; x < dim; x++) std::cout<<xiNode[x]<<"\t";
//                std::cout<<std::endl;
//#endif
                for (auto x = 0; x < 3; x++) xNode[x] = nodeCoords.GetVal(iNode, x);
                gmesh->NodeVec()[newIndex].Initialize(xNode, *gmesh);
            }
            TPZVec<int64_t> nodes(nNodes, 0);
            for (auto i = 0; i < nNodes; i++) nodes[i] = i;
            int64_t index;
            TPZGeoEl *gel = gmesh->CreateGeoElement(elType, nodes, 1, index, 0);

            //computing directions
            TPZFMatrix<REAL> gradX(3, dim, 0);
            for (auto x = 0; x < dim; x++) gradX(x, x) = 1;
            const int nVec = dim * nSides;
            TPZFMatrix<REAL> masterDirections(3, nVec), deformedDirections(3, nVec);
            TTopol::ComputeHCurlDirections(gradX, masterDirections);//these are the vectors on the master element
            TPZVec<REAL> xiCenter(dim, 0);
            TPZFMatrix<REAL> jacobian(dim, dim, 0), axes(dim, 3), jacInv(dim, dim, 0);
            REAL detJac = 0;
            gel->JacobianXYZ(xiCenter, jacobian, axes, detJac, jacInv);

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
            //calculating tangent vectors for all edges
            TPZFMatrix<REAL> edgeTangentVecs(nEdges,3);
            TPZVec<REAL> edgeLengthVec(nEdges,0);
            for (auto iEdge = 0; iEdge < nEdges; iEdge++) {
                const int edgeIndex = nNodes + iEdge;

                TPZManVector<int, 2> edgeNodes(2, 0);
                for (auto i = 0; i < 2; i++) edgeNodes[i] = gel->SideNodeIndex(edgeIndex, i);

                TPZVec<REAL> edgeTgVector(3, 0);
                REAL edgeLength = 0;
                {
                    TPZVec<REAL> p0(3), p1(3);
                    gmesh->NodeVec()[edgeNodes[0]].GetCoordinates(p0);
                    gmesh->NodeVec()[edgeNodes[1]].GetCoordinates(p1);
                    edgeLength =
                            std::sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                                      (p0[2] - p1[2]) * (p0[2] - p1[2]));
                    for (auto x = 0; x < 3; x++) edgeTangentVecs(iEdge,x) = (p1[x] - p0[x]) / edgeLength;
                    edgeLengthVec[iEdge] = edgeLength;
                }
            }
            //testing directions associated with edges
            for (auto iEdge = 0; iEdge < nEdges; iEdge++) {
                std::cout << "\tedge " << iEdge+1<< " out of " << nEdges << std::endl;
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
                std::cout << "\tedgeIndex " << edgeIndex <<"\ttesting vea vectors" << std::endl;
                for(auto iVec = 0; iVec <  2; iVec++) {// 2 v^{e,a}
                    std::cout<<"\t\tvector  "<<iVec<<" out of "<<2<<std::endl;
                    trace = 0;
                    for(auto x = 0; x < 3; x++) trace += edgeTgVector[x] * deformedDirections(x,2*iEdge+iVec);
                    trace *=edgeLength;
    #ifdef NOISY_HCURL
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

            //auxiliary lambda expressions
            auto VectorialProduct = [](TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result){
                REAL x1=v1[0], y1=v1[1],z1=v1[2];
                REAL x2=v2[0], y2=v2[1],z2=v2[2];
                result.Resize(v1.NElements());
                result[0]=y1*z2-z1*y2;
                result[1]=z1*x2-x1*z2;
                result[2]=x1*y2-y1*x2;
            };

            auto ComputeNormal = [VectorialProduct,gmesh] (TPZVec<int> faceNodes,TPZVec<REAL> &result){
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
                VectorialProduct(v1,v2,result);
                REAL norm = 0;
                for(auto x = 0; x < 3; x++) norm += result[x] * result[x];
                norm = std::sqrt(norm);
                for(auto x = 0; x < 3; x++) result[x] /= norm;
            };
            //now the testing for directions associated with faces
            //@TODOFran: how to check the trace of face vectors?
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

            for(auto iFace = 0; iFace < nFaces; iFace++){
                std::cout<<"\tface "<<iFace+1<<" out of "<<nFaces<<std::endl;
                const int faceIndex = nNodes + nEdges + iFace;
                TPZManVector<int,3> faceNodes(3,0);
                for(auto i = 0; i < 3; i++) faceNodes[i] = gel->SideNodeIndex(faceIndex,i);
                TPZManVector<REAL,3> deformedFaceNormal(3,0);
                ComputeNormal(faceNodes,deformedFaceNormal);
                const REAL faceArea = gel->SideArea(faceIndex);

                std::cout<<"\t\ttesting vfe vectors"<<std::endl;
                const int nFaceEdges = faceEdges[iFace].size();
                for(auto iVec = 0; iVec < nFaceEdges; iVec++) {
                    std::cout<<"\t\t\tvec "<<iVec<<" out of "<<nFaceEdges<<std::endl;

                    TPZManVector<REAL,3> vfe(3,0);
                    const int vfeIndex = firstVfeVec[iFace]+iVec;
                    for(auto x = 0; x < 3; x++) vfe[x] = deformedDirections(x,vfeIndex);
                    TPZManVector<REAL,3> temp(3,0);
                    TPZManVector<REAL,3> tangentialTrace(3,0);
                    VectorialProduct(deformedFaceNormal,vfe,temp);
                    VectorialProduct(deformedFaceNormal,temp,tangentialTrace);
#ifdef NOISY_HCURL
                    std::cout<<"\t\t\tdeformed vfe:"<<std::endl;
                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< vfe[x]<<"\t";
                    std::cout<<std::endl;
                    std::cout<<"\t\t\ttangential trace over edge:"<<std::endl<<"\t\t";
                    for(auto x = 0; x < 3; x++) std::cout<<tangentialTrace[x]<<"\t";
                    std::cout<<std::endl;
#endif
//                    REAL traceNorm = 0;
//                    for(auto x = 0; x < 3; x++) traceNorm +=  tangentialTrace[x] * tangentialTrace[x];
//                    traceNorm = sqrt(traceNorm) * faceArea;
//                    bool checkTrace =
//                            std::abs(traceNorm - 1) < tol;
//                    BOOST_CHECK(checkTrace);
                    if(nFaces == 1) continue;
                    //check that this vector is normal to the edge \hat{e} adjacent by a to e
                    const int edgeIndex = faceEdges[iFace][iVec];
                    int neighFaceIndex = -1;
                    for(int iNeigh = 0; iNeigh < nFaces; iNeigh++){
                        const int candidateFaceIndex = iNeigh + nNodes + nEdges;
                        const int nNeighFaceEdges = faceEdges[iNeigh].size();
                        for(int iEdge = 0; iEdge < nNeighFaceEdges; iEdge++){
                            const bool cond1 = faceEdges[iNeigh][iEdge] == edgeIndex;
                            const bool cond2 = candidateFaceIndex != faceIndex;
                            if(cond1 && cond2){
                                neighFaceIndex = candidateFaceIndex;
                            }
                        }
                    }
                    if(neighFaceIndex == -1){
                        DebugStop();
                    }
                    else{
                        TPZManVector<int,3> neighFaceNodes(3,0);
                        for(auto i = 0; i < 3; i++) neighFaceNodes[i] = gel->SideNodeIndex(neighFaceIndex,i);
                        TPZManVector<REAL,3> neighNormal(3,0);
                        ComputeNormal(neighFaceNodes,neighNormal);
                        TPZManVector<REAL,3> temp(3,0);
                        TPZManVector<REAL,3> neighTangentialTrace(3,0);
                        VectorialProduct(neighNormal,vfe,temp);
                        VectorialProduct(neighNormal,temp,neighTangentialTrace);
                        REAL neighTraceNorm = 0;
                        for(auto x = 0; x < 3; x++) neighTraceNorm +=  neighTangentialTrace[x] * neighTangentialTrace[x];
                        bool checkNeighTrace =
                                std::abs(neighTraceNorm) < tol;

                        if(!checkNeighTrace){
                            std::cout<<"\t\t\tERROR on direction "<<vfeIndex<<std::endl;
                            std::cout<<"\t\t\ttangential trace over neighbour face: "<<std::endl<<"\t\t";
                            for (auto x = 0; x < 3; x++) std::cout << neighTangentialTrace[x]<< "\t";
                            std::cout<<std::endl;
                            std::cout<<"\t\t\tneighbour face: "<<neighFaceIndex<<std::endl;
                            std::cout<<"\t\t\tneighbour face normal (deformed): "<<std::endl<<"\t\t";
                            for (auto x = 0; x < 3; x++) std::cout << neighNormal[x]<< "\t";
                            std::cout<<std::endl;
                        }
                        BOOST_CHECK(checkNeighTrace);
                    }

                }
//                std::cout<<"\t\ttesting vft vectors"<<std::endl;
//                for(auto iVec = 0; iVec < 2; iVec++) {//there are two v^{F,T} vectors per face
//                    std::cout<<"\tdirection "<<iVec<<" out of 2"<<std::endl;
//                    TPZManVector<REAL,3> vft(3,0);
//                    for(auto x = 0; x < 3; x++) vft[x] = deformedDirections(x,firstVftVec + iFace * 2 + iVec);
//                    TPZManVector<REAL,3> temp(3,0);
//                    TPZManVector<REAL,3> tangentialTrace(3,0);
//                    VectorialProduct(deformedFaceNormal,vft,temp);
//                    VectorialProduct(deformedFaceNormal,temp,tangentialTrace);
//
//#ifdef NOISY_HCURL
//                    std::cout<<"\t\tdeformed vft:"<<std::endl;
//                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< vft[x]<<"\t";
//                    std::cout<<std::endl;
//                    std::cout<<"\t\ttangential trace over edge:"<<std::endl;
//                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< tangentialTrace[x]<<"\t";
//                    std::cout<<std::endl;
//#endif
//                    REAL traceNorm = 0;
//                    for(auto x = 0; x < 3; x++) traceNorm +=  tangentialTrace[x] * tangentialTrace[x];
//                    traceNorm = sqrt(traceNorm) * faceArea;
//                    bool checkTrace =
//                            (traceNorm - 1) < tol;
//                    BOOST_CHECK(checkTrace);
//                }
            }
        }//hcurltest::CompareTraces




    }//namespace


BOOST_AUTO_TEST_SUITE_END()

#ifdef NOISY_HCURL
#undef NOISY_HCURL
#endif
#ifdef NOISYVTK_HCURL
#undef NOISYVTK_HCURL
#endif

#endif

