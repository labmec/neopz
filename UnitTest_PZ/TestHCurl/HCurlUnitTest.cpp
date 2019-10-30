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

//#define NOISY_HCURL //outputs useful debug info
//#define NOISYVTK_HCURL//outputs even more debug info, printing relevant info in .vtk format

std::string dirname = PZSOURCEDIR;

#ifdef USING_BOOST


BOOST_AUTO_TEST_SUITE(hcurl_tests)
    namespace hcurltest{
        constexpr REAL tol = 1e-10;
        template <class TTopol>
        void ComparePermutedVectors();
    }

    BOOST_AUTO_TEST_CASE(hcurl_permutation_tests) {
        InitializePZLOG();
        hcurltest::ComparePermutedVectors<pztopology::TPZTriangle>();
        hcurltest::ComparePermutedVectors<pztopology::TPZQuadrilateral>();
    }

    namespace hcurltest{


        template <class TTopol>
        void ComparePermutedVectors() {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            auto elType = TTopol::Type();

            const int64_t nSides = TTopol::NumSides();
            const int64_t nNodes = TTopol::NumSides(0);
            const int64_t nEdges = TTopol::NumSides(1);
            const int64_t nFaces = TTopol::NumSides(2);
            const int64_t dim = TTopol::Dimension;

            TPZGeoMesh *gmesh = new TPZGeoMesh();
            ///creating the mesh nodes. they will coincide with the master element nodes
            TPZVec<REAL> nodeLocalCoord(0), xiNode(dim,0);
            int64_t newIndex = -1;
            for(int iNode = 0; iNode < nNodes; iNode++){
                auto transf = TTopol::TransformSideToElement(iNode);
                transf.Apply(nodeLocalCoord,xiNode);
                newIndex = gmesh->NodeVec().AllocateNewElement();
//#ifdef NOISY_HCURL
//                std::cout<<"node "<<iNode<<" coords:"<<std::endl;
//                for(auto x = 0; x < dim; x++) std::cout<<xiNode[x]<<"\t";
//                std::cout<<std::endl;
//#endif
                gmesh->NodeVec()[newIndex].Initialize(xiNode, *gmesh);
            }

            const int64_t nPermutations = TTopol::NPermutations;

            auto CreatePermutedEl = [gmesh,elType,nNodes]( TPZVec<int> &currentPermutation){
                TPZVec<int64_t> nodesPerm(nNodes,0);
                for(auto x = 0 ; x < nNodes; x++) nodesPerm[x] = currentPermutation[x];
                int64_t index, matId = 1;
                return gmesh->CreateGeoElement(elType,nodesPerm,matId,index,0);
            };

            int permute = 0;
            TPZVec<int> currentPermutation(nSides,0);
            pztopology::GetPermutation<TTopol>(permute, currentPermutation);
            TPZGeoEl * originalEl = CreatePermutedEl(currentPermutation);

            TPZFMatrix<REAL> gradX(3,dim,0);
            for(auto x = 0 ; x < dim; x++) gradX(x,x) = 1;
            const int nVec = dim * nSides;
            TPZFMatrix<REAL> directions(3,nVec);
            TTopol::ComputeHCurlDirections(gradX,directions);//these are the vectors on the master element



            TPZGeoEl *permutedEl = nullptr;
            TPZFMatrix<REAL> permutedJacobian(dim,dim,0),axes(dim,3), permutedJacInv(dim,dim,0);
            REAL detPermutedJac = 0;
            TPZVec<REAL> xiCenter(dim,0);
            TTopol::CenterPoint(nSides - 1, xiCenter);
            TPZFMatrix<REAL> permutedDirections(3,nVec,0);
            for(permute = 1; permute < nPermutations; permute++){
                pztopology::GetPermutation<TTopol>(permute, currentPermutation);
                permutedEl = CreatePermutedEl(currentPermutation);
                //in all valid permutations the gradX is constant therefore the xi point should not matter
                permutedEl->JacobianXYZ(xiCenter,permutedJacobian,axes,detPermutedJac,permutedJacInv);

                #ifdef NOISY_HCURL
                std::cout<<"permutation "<<permute<<" out of "<<nPermutations<<". side ordering:"<<std::endl;
                for(auto iSide = 0; iSide < nSides; iSide++) std::cout<<currentPermutation[iSide]<<"\t";
                std::cout<<std::endl;
                permutedJacobian.Print("\tJacobian:");
                axes.Print("\tAxes:");
                #endif
                for(auto iVec = 0; iVec <  nVec; iVec++){
                    //covariant piola transform: J^{-T}
                    TPZManVector<REAL,3> tempDirection(dim,0);
                    for(auto i = 0; i < dim; i++){
                        tempDirection[i] = 0;
                        for(auto j = 0; j< dim; j++)    tempDirection[i] += permutedJacInv(j,i) * directions(j,iVec);
                    }
                    for(auto i = 0; i < 2; i++){
                        permutedDirections(i,iVec) = 0;
                        for(auto j = 0; j< dim; j++)    permutedDirections(i,iVec) += axes(j,i) * tempDirection[j];
                    }
                }
                //testing directions associated with edges
                for(auto iEdge = 0; iEdge < nEdges; iEdge++){
                    const int originalEdgeIndex = nNodes + iEdge;
                    int permutedIEdge = -1;
                    for(auto pEdge = 0; pEdge < nEdges; pEdge++){
                        if(currentPermutation[pEdge+nNodes] == originalEdgeIndex){
                            permutedIEdge = pEdge;
                            break;
                        }
                    };
                    const int permutedEdgeIndex = permutedIEdge + nNodes;
                    TPZVec<REAL> edgeTgVector(3,0);
                    REAL tgNorm = 0;
                    for(auto x = 0; x < 3; x++) {
                        edgeTgVector[x] = directions(x,2*nEdges+iEdge);
                        tgNorm += edgeTgVector[x] * edgeTgVector[x];
                    }
                    tgNorm = std::sqrt(tgNorm);
                    for(auto x = 0; x < 3; x++) edgeTgVector[x] /= tgNorm;
                    REAL sideOrient = 0;
                    for(auto x = 0; x < 3; x++) sideOrient += edgeTgVector[x]* permutedDirections(x,2*nEdges+permutedIEdge);
                    sideOrient /= std::abs(sideOrient);

                    REAL diff = 0;
                    for(auto x = 0; x < 3; x++) diff += directions(x,2*nEdges+iEdge) - sideOrient * permutedDirections(x,2*nEdges+permutedIEdge);
                    bool areTangentialVectorsDifferent =
                            std::abs(diff) < tol;
                    BOOST_CHECK(areTangentialVectorsDifferent);
                    #ifdef NOISY_HCURL
                    std::cout<<"\tedge "<<originalEdgeIndex<<std::endl;
                    std::cout<<"\tpermuted edge "<<permutedEdgeIndex<<std::endl;
                    std::cout<<"\tside orient "<<sideOrient<<std::endl;
                    std::cout<<"\ttesting ve vectors"<<std::endl;
                    std::cout<<"\t\toriginal:"<<std::endl;
                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< directions(x,2*nEdges+iEdge)<<"\t";
                    std::cout<<"\n\t\tpermuted:"<<std::endl;
                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< sideOrient * permutedDirections(x,2*nEdges+permutedIEdge)<<"\t";
                    std::cout<<std::endl;

                    std::cout<<"\t\ttangent vector: ";
                    for(auto x = 0; x < 3; x++) std::cout<<edgeTgVector[x]<<"\t";
                    std::cout<<std::endl;
                    std::cout<<"\t\ttesting vea vectors"<<std::endl;
                    #endif
                    REAL originalTangentialTrace = 0, permutedTangentialTrace = 0;
                    for(auto iVec = 0; iVec <  2; iVec++) {// 2 v^{e,a}
                        originalTangentialTrace = 0;
                        for(auto x = 0; x < 3; x++) originalTangentialTrace += edgeTgVector[x] * directions(x,2*iEdge+iVec);
                        permutedTangentialTrace = 0;
                        for(auto x = 0; x < 3; x++) permutedTangentialTrace += edgeTgVector[x] * sideOrient * permutedDirections(x,2*permutedIEdge+iVec);
                        #ifdef NOISY_HCURL
                        std::cout<<"\tdirection "<<iVec<<std::endl;
                        std::cout<<"\t\toriginal:"<<std::endl;
                        for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< directions(x,2*iEdge+iVec)<<"\t";
                        std::cout<<std::endl;
                        std::cout<<"\t\ttangential trace over edge: "<<originalTangentialTrace<<std::endl;
                        std::cout<<"\t\tpermuted:"<<std::endl;
                        for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< sideOrient * permutedDirections(x,2*permutedIEdge+iVec)<<"\t";
                        std::cout<<std::endl;
                        std::cout<<"\t\ttangential trace over edge: "<<permutedTangentialTrace<<std::endl;
                        #endif

                        bool areTracesDifferent =
                                std::abs(originalTangentialTrace - permutedTangentialTrace) < tol;
                        BOOST_CHECK(areTracesDifferent);
                    }
                }
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
                    TPZVec<REAL> p1(3),p2(3),p3(3);
                    gmesh->NodeVec()[faceNodes[0]].GetCoordinates(p1);
                    gmesh->NodeVec()[faceNodes[1]].GetCoordinates(p2);
                    gmesh->NodeVec()[faceNodes[2]].GetCoordinates(p3);
                    v1[0]=p1[0]-p2[0];
                    v1[1]=p1[1]-p2[1];
                    v1[2]=p1[2]-p2[2];
                    v2[0]=p3[0]-p2[0];
                    v2[1]=p3[1]-p2[1];
                    v2[2]=p3[2]-p2[2];
                    VectorialProduct(v1,v2,result);
                    REAL norm = 0;
                    for(auto x = 0; x < 3; x++) norm += result[x] * result[x];
                    norm = std::sqrt(norm);
                    for(auto x = 0; x < 3; x++) result[x] /= norm;
                };
                const int nEdgeVectors = nEdges * 3;
                //testing directions associated with faces
                for(auto iFace = 0; iFace < nFaces; iFace++){
                    const int originalFaceIndex = nNodes + nEdges + iFace;
                    int permutedIFace = -1;
                    for(auto pFace = 0; pFace < nEdges; pFace++){
                        if(currentPermutation[pFace+nEdges+nNodes] == originalFaceIndex){
                            permutedIFace = pFace;
                            break;
                        }
                    };
                    const int permutedFaceIndex = permutedIFace + nNodes + nEdges;

                    TPZManVector<int,3> faceNodes(3,0);
                    for(auto i = 0; i < 3; i++) faceNodes[i] = originalEl->SideNodeIndex(originalFaceIndex,i);
                    TPZVec<REAL> originalFaceNormalVector(3,0);
                    ComputeNormal(faceNodes,originalFaceNormalVector);

                    for(auto i = 0; i < 3; i++) faceNodes[i] = permutedEl->SideNodeIndex(permutedFaceIndex,i);
                    TPZVec<REAL> permutedFaceNormalVector(3,0);
                    ComputeNormal(faceNodes,permutedFaceNormalVector);

                    REAL sideOrient = 0;
                    for(auto x = 0; x < 3; x++) sideOrient += originalFaceNormalVector[x] * permutedFaceNormalVector[x];
                    sideOrient /= std::abs(sideOrient);

//                    #ifdef NOISY_HCURL
//                    std::cout<<"\tface "<<originalFaceIndex<<std::endl;
//                    std::cout<<"\tpermuted face "<<permutedFaceIndex<<std::endl;
//                    std::cout<<"\tside orient "<<sideOrient<<std::endl;
//                    std::cout<<"\ttesting ve vectors"<<std::endl;
//                    std::cout<<"\t\toriginal:"<<std::endl;
//                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< directions(x,2*nEdges+iEdge)<<"\t";
//                    std::cout<<"\n\t\tpermuted:"<<std::endl;
//                    for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< sideOrient * permutedDirections(x,2*nEdges+permutedIEdge)<<"\t";
//                    std::cout<<std::endl;
//
//                    std::cout<<"\t\ttangent vector: ";
//                    for(auto x = 0; x < 3; x++) std::cout<<edgeTgVector[x]<<"\t";
//                    std::cout<<std::endl;
//                    std::cout<<"\t\ttesting vea vectors"<<std::endl;
//                    #endif
//                    REAL originalTangentialTrace = 0, permutedTangentialTrace = 0;
//                    for(auto iVec = 0; iVec <  2; iVec++) {// 2 v^{e,a}
//                        originalTangentialTrace = 0;
//                        for(auto x = 0; x < 3; x++) originalTangentialTrace += edgeTgVector[x] * directions(x,2*iEdge+iVec);
//                        permutedTangentialTrace = 0;
//                        for(auto x = 0; x < 3; x++) permutedTangentialTrace += edgeTgVector[x] * sideOrient * permutedDirections(x,2*permutedIEdge+iVec);
//                        #ifdef NOISY_HCURL
//                        std::cout<<"\tdirection "<<iVec<<std::endl;
//                        std::cout<<"\t\toriginal:"<<std::endl;
//                        for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< directions(x,2*iEdge+iVec)<<"\t";
//                        std::cout<<std::endl;
//                        std::cout<<"\t\ttangential trace over edge: "<<originalTangentialTrace<<std::endl;
//                        std::cout<<"\t\tpermuted:"<<std::endl;
//                        for(auto x = 0; x < 3; x++) std::cout<<"\t\t"<< sideOrient * permutedDirections(x,2*permutedIEdge+iVec)<<"\t";
//                        std::cout<<std::endl;
//                        std::cout<<"\t\ttangential trace over edge: "<<permutedTangentialTrace<<std::endl;
//                        #endif
//
//                        bool areTracesDifferent =
//                                std::abs(originalTangentialTrace - permutedTangentialTrace) < tol;
//                        BOOST_CHECK(areTracesDifferent);
//                    }
                }
            }
            delete gmesh;
        }



    }


BOOST_AUTO_TEST_SUITE_END()

#ifdef NOISY_HCURL
#undef NOISY_HCURL
#endif
#ifdef NOISYVTK_HCURL
#undef NOISYVTK_HCURL
#endif

#endif

