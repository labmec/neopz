//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzmanvector.h"
#include "pztrnsform.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"
#include "pzgeotetrahedra.h"
#include "pzgeotriangle.h"

#include <catch2/catch.hpp>

using namespace pztopology;
namespace topologytests{
    template <class top>
    void TestingSideProjections();
    template <class top>
    void TestingSideNodeProjections();
    template <class top>
    void TestingConstantDivergent();
}

TEST_CASE("projection_tests_1","[topology_tests]")
{
    topologytests::TestingSideProjections<TPZLine>();
    topologytests::TestingSideProjections<TPZTriangle>();
    topologytests::TestingSideProjections<TPZQuadrilateral>();
    topologytests::TestingSideProjections<TPZTetrahedron>();
    topologytests::TestingSideProjections<TPZCube>();
    topologytests::TestingSideProjections<TPZPrism>();
    topologytests::TestingSideProjections<TPZPyramid>();
}

TEST_CASE("projection_tests_2","[topology_tests]")
{
    topologytests::TestingSideNodeProjections<TPZTriangle>();
    topologytests::TestingSideNodeProjections<TPZQuadrilateral>();
    topologytests::TestingSideNodeProjections<TPZTetrahedron>();
    topologytests::TestingSideNodeProjections<TPZCube>();
    topologytests::TestingSideNodeProjections<TPZPrism>();
    topologytests::TestingSideNodeProjections<TPZPyramid>();
}

TEST_CASE("constant_divergent_test","[topology_tests]")
{
    // topologytests::TestingConstantDivergent<TPZTriangle>();
    topologytests::TestingConstantDivergent<TPZQuadrilateral>();
    // topologytests::TestingConstantDivergent<TPZTetrahedron>();
    // topologytests::TestingConstantDivergent<TPZCube>();
    // topologytests::TestingConstantDivergent<TPZPrism>();
}

namespace topologytests{
    template <class top>
    void TestingSideProjections() {
        static std::string testName = __PRETTY_FUNCTION__;
        auto type = top::Type();
        int dim = top::Dimension;
        int nsides = top::NSides;
        for (auto iSide=0; iSide<nsides; iSide++) {
            TPZTransform<> tr1 = top::TransformSideToElement(iSide);
            TPZTransform<> tr2 = top::TransformElementToSide(iSide);
            TPZTransform<> tr3;
            tr3 = tr2.Multiply(tr1);
            TPZMatrix<REAL>& mat = tr3.Mult();
            const REAL tol = GetTolerance();
            const int nRows = mat.Rows();
            const int nCols = mat.Cols();
            for(int i = 0; i < nRows; i++){
                for(int j = 0; j < nCols; j++){
                    bool cond;
                    if (j == i) {
                        cond = mat(i,j)-1 < tol;
                    }else{
                        cond = std::abs(mat(i,j) ) < tol;
                    }
                    if(!cond){
                      std::cerr
                          << "\n" + testName + " failed" +
                                 "\ntopology: " + MElementType_Name(type) +
                                 "\n" + "side: " + std::to_string(iSide);
                    }
                    REQUIRE(cond);
                }
            }
        }
    }//TestingSideProjections

    template <class top>
    void TestingSideNodeProjections() {
        static std::string testName = __PRETTY_FUNCTION__;
        auto type = top::Type();
        int dim = top::Dimension;
        int nSides = top::NSides;
        int nNodes = top::NCornerNodes;
        TPZFNMatrix<24,REAL> cornerNodes(nNodes,dim);
        TPZManVector<REAL,3> node(dim);
        for (auto iNode=0; iNode<nNodes; iNode++) {
            top::ParametricDomainNodeCoord(iNode, node);
            for(auto x = 0; x < dim; x++) cornerNodes(iNode,x) = node[x];
        }
        const REAL tol = GetTolerance();
        for (auto iSide=nNodes; iSide<nSides - 1; iSide++) {
            const auto sideDim = top::SideDimension(iSide);
            TPZManVector<REAL,3> paramNode(sideDim);
            TPZTransform<> tr = top::TransformSideToElement(iSide);
            const auto nSideNodes = top::NSideNodes(iSide);
            for(auto iNode = 0; iNode < nSideNodes; iNode++){
                switch(nSideNodes){
                    case 2://edge
                        TPZLine::ParametricDomainNodeCoord(iNode,paramNode);
                        break;
                    case 3://triangle
                        TPZTriangle::ParametricDomainNodeCoord(iNode,paramNode);
                        break;
                    case 4://quadrilateral
                        TPZQuadrilateral::ParametricDomainNodeCoord(iNode,paramNode);
                        break;
                    default:
                        DebugStop();
                }
                tr.Apply(paramNode,node);
                const auto nodeId = top::SideNodeLocId(iSide,iNode);
                REAL diff = 0;
                for(auto x = 0; x < sideDim; x++) diff += (node[x] - cornerNodes(nodeId,x)) * (node[x] - cornerNodes(nodeId,x));
                bool cond = std::sqrt(diff) < tol;
                if(!cond){
                  std::cerr << "\n" + testName + " failed" +
                                   "\ntopology: " + MElementType_Name(type) +
                                   "\n" + "side: " + std::to_string(iSide);
                }
                REQUIRE(cond);
            }
        }
    }//TestingSideNodeProjections


    template <class top>
    void TestingConstantDivergent() {
        
        static std::string testName = __PRETTY_FUNCTION__;
        const  int64_t pOrderIntRule = 6;
        const auto nSides = top::NSides;
        const auto nCorner = top::NCornerNodes;
        const auto dim = top::Dimension;
        const auto nFaces = top::NFacets;
        auto type = top::Type();
        const REAL tol = 1.e-6;

        auto gel = pzgeom::TPZNodeRep<nCorner,top>();

        TPZIntPoints* intRule = gel.CreateSideIntegrationRule(nSides-1, pOrderIntRule);     

        TPZFMatrix<REAL> vecDiv(dim,nFaces);
        TPZFMatrix<REAL> vecCurl(dim,nFaces);
        TPZFMatrix<REAL> RT0Function(dim,nFaces);
        TPZVec<REAL> div(nFaces);
        TPZVec<REAL> node(dim);
        
        //For HCurl
        TPZVec<int> transformIds(nSides-nCorner,1);
        TPZVec<int64_t> nodes(nCorner, 0);
        for (auto i = 0; i < nCorner; i++) nodes[i] = gel.fNodeIndexes[i];
        for(auto iSide = 0 ; iSide < nSides - nCorner; iSide++){
            transformIds[iSide] = top::GetTransformId(nCorner + iSide, nodes);
        }

        const int npts = intRule->NPoints();
        REAL w;
        for (auto ipt = 0; ipt < npts; ipt++) {
            intRule->Point(ipt, node, w);

            TPZFNMatrix<nCorner> phis(nCorner,1);
            TPZFNMatrix<nSides*dim*dim> directionsHDiv(3,nSides*dim);
            TPZFNMatrix<nSides*dim*dim> directionsHCurl(3,nSides*dim);
            
            TPZFNMatrix<dim*dim> gradx(3,3);
            TPZFNMatrix<dim*nCorner> dphis(dim,nCorner);
            gradx.Identity();
            vecDiv.Zero();
            vecCurl.Zero();
            RT0Function.Zero();
            div.Fill(0.);

            top::ComputeHDivDirections(gradx,directionsHDiv);
            top::ComputeHCurlDirections(gradx,directionsHCurl,transformIds);
            top::Shape(node,phis,dphis);

            int first_face = nSides-1-nFaces;
            int firstVecIndex = 0;
            for (size_t iface = first_face; iface < nSides-1; iface++)
            {
                int face_count = iface - first_face;
                int nsubsides = top::NContainedSides(iface);
                int ncorner = top::NSideNodes(iface);

                for (size_t ivec = 0; ivec < ncorner; ivec++)
                {
                    TPZManVector<REAL,dim> vec(dim);
                    int vecIndex = firstVecIndex + ivec;
                    int vertex = top::SideNodeLocId(iface,ivec);
                    REAL divlocal = 0.;
                    REAL divlocal2 = 0.;

                    for (size_t i = 0; i < dim; i++)
                    {
                        divlocal += directionsHDiv(i,vecIndex) * dphis(i,vertex) / top::NSideNodes(iface);
                        vecDiv(i,face_count) += directionsHDiv(i,vecIndex) / top::NSideNodes(iface);
                        vecCurl(i,face_count) += directionsHCurl(i,vecIndex) / top::NSideNodes(iface);
                    }//i
                    div[face_count] += divlocal;
                }
                firstVecIndex += nsubsides;
            }

            //Compute the divergent for each face
            top::ComputeConstantHDiv(node,RT0Function,div);

            std::cout << "Point = " << node << std::endl;
            std::cout << "Vecdiv = " << vecCurl << std::endl;
            // std::cout << "RT0 = " << RT0Function << std::endl;
            // std::cout << "DIV = " << div << std::endl;
            
            
            // // Checks if all faces have the same divergent value.
            // bool cond = true;
            // for (int i = 0; i < dim; i++){
            //     for (int j = 0; j < nFaces; j++){
            //         if (fabs(vecDiv(i,j)-RT0Function(i,j)) > tol){
            //             cond = false;
            //         } 
            //     }
            // }
            // if(!cond){
            //     std::cerr << "\n" + testName + " failed" +
            //                     "\ntopology: " + MElementType_Name(type);
            // }
            // REQUIRE(cond);
        }
    }//Testing Constant Divergent

}//namespace