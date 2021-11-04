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
#include "TPZShapeHCurl.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeData.h"
#include "TPZMaterialData.h"
#include "pzelctemp.h"
#include "TPZShapeHDivKernel.h"
#include "TPZShapeHDivKernel2D.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapetetra.h"



#include <catch2/catch.hpp>

using namespace pztopology;
using namespace pzshape;
namespace topologytests{
    template <class top>
    void TestingSideProjections();
    template <class top>
    void TestingSideNodeProjections();
    template <class top>
    void TestingConstantDivergent();
    template <class top,class TSHAPE>
    void TestingConstantCurl();
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
    topologytests::TestingConstantDivergent<TPZTriangle>();
    topologytests::TestingConstantDivergent<TPZQuadrilateral>();
    topologytests::TestingConstantDivergent<TPZTetrahedron>();
    topologytests::TestingConstantDivergent<TPZCube>();
    topologytests::TestingConstantDivergent<TPZPrism>();
}

TEST_CASE("constant_curl_test_new","[topology_tests]")
{
    topologytests::TestingConstantCurl<TPZTriangle,TPZShapeTriang>();
    topologytests::TestingConstantCurl<TPZQuadrilateral,TPZShapeQuad>();
    topologytests::TestingConstantCurl<TPZTetrahedron,TPZShapeTetra>();
    topologytests::TestingConstantCurl<TPZCube,TPZShapeCube>();
    topologytests::TestingConstantCurl<TPZPrism,TPZShapePrism>();
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
        const  int64_t pOrderIntRule = 5;
        const auto nSides = top::NSides;
        const auto nCorner = top::NCornerNodes;
        const auto dim = top::Dimension;
        const auto nFaces = top::NFacets;
        auto type = top::Type();
        const REAL tol = 1.e-6;

        auto gel = pzgeom::TPZNodeRep<nCorner,top>();

        TPZIntPoints* intRule = gel.CreateSideIntegrationRule(nSides-1, pOrderIntRule);     

        TPZFMatrix<REAL> vecDiv(dim,nFaces);
        TPZFMatrix<REAL> RT0Function(dim,nFaces);
        TPZVec<REAL> div(nFaces);
        TPZVec<REAL> divRT0(nFaces);
        TPZVec<REAL> node(dim);

        const int npts = intRule->NPoints();
        REAL w;

        for (auto ipt = 0; ipt < npts; ipt++) {
            intRule->Point(ipt, node, w);

            TPZFNMatrix<nCorner> phis(nCorner,1);
            TPZFNMatrix<nSides*dim*dim> directionsHDiv(3,nSides*dim);            
            TPZFNMatrix<dim*dim> gradx(3,3);
            TPZFNMatrix<dim*nCorner> dphis(dim,nCorner);
            gradx.Identity();
            vecDiv.Zero();
            RT0Function.Zero();
            div.Fill(0.);

            top::ComputeHDivDirections(gradx,directionsHDiv);
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

                    for (size_t i = 0; i < dim; i++)
                    {
                        divlocal += directionsHDiv(i,vecIndex) * dphis(i,vertex) / top::NSideNodes(iface);
                        vecDiv(i,face_count) += directionsHDiv(i,vecIndex) * phis(vertex) / top::NSideNodes(iface);
                    }//i
                    div[face_count] += divlocal;
                }
                firstVecIndex += nsubsides;
            }

            divRT0.Fill(0.);
            //Compute the divergent for each face
            top::ComputeConstantHDiv(node,RT0Function,divRT0);  
            
            // Checks if all faces have the same divergent value.
            bool condHdiv = true;
            bool condHdivRT0 = true;
            for (int i = 0; i < dim; i++){
                for (int j = 0; j < nFaces; j++){
                    if (fabs(vecDiv(i,j)-RT0Function(i,j)) > tol){
                        condHdiv = false;
                    } 
                }
            }
            for (int j = 0; j < nFaces; j++){
                if (fabs(div[j]-divRT0[j]) > tol){
                    condHdivRT0 = false;
                } 
            }
            if(!condHdiv || !condHdivRT0){
                std::cerr << "\n" + testName + " failed Hdiv " +
                                "\ntopology: " + MElementType_Name(type);
            }
            REQUIRE(condHdiv);
            REQUIRE(condHdivRT0);
        }
    }//Testing Constant Divergent



    template <class top,class TSHAPE>
    void TestingConstantCurl() {
        
        static std::string testName = __PRETTY_FUNCTION__;
        const  int64_t pOrderIntRule = 5;
        const auto nSides = top::NSides;
        const auto nCorner = top::NCornerNodes;
        const auto dim = top::Dimension;
        const auto nFaces = top::NFacets;
        auto type = top::Type();
        const REAL tol = 1.e-6;
        auto gel = pzgeom::TPZNodeRep<nCorner,top>();

        TPZIntPoints* intRule = gel.CreateSideIntegrationRule(nSides-1, pOrderIntRule);     
        int nEdges = 0;
        if (dim == 2){
            nEdges = nFaces;
        } else if (dim == 3){
            nEdges = nSides - 1 - nCorner - nFaces;
        }
        
        TPZVec<REAL> node(dim);
                
        TPZShapeData shapedata;
        TPZManVector<int64_t,nCorner> ids(nCorner,0);
        for(auto i=0; i<nCorner; i++) ids[i] = i;        
        auto &conOrders = shapedata.fHDivConnectOrders;
        constexpr auto nConnects = nSides - nCorner;
        conOrders.Resize(nConnects,-1);
        for(auto i = 0; i < nConnects; i++) conOrders[i] = 1;

        TPZShapeHCurl<TSHAPE>::Initialize(ids, conOrders, shapedata);

        const int npts = intRule->NPoints();
        auto nshape = shapedata.fSDVecShapeIndex.size();

        for (auto ipt = 0; ipt < npts; ipt++) {
            REAL w;
            intRule->Point(ipt, node, w);

            TPZFMatrix<REAL> phiHCurl(3,nshape);
            TPZFMatrix<REAL> curlHCurl(3,nshape);
            TPZFMatrix<REAL> N0Function(dim,nEdges);
            TPZFMatrix<REAL> curlN0(3,nEdges);
            N0Function.Zero();
            curlN0.Zero();
            phiHCurl.Zero();
            curlHCurl.Zero();

            TPZShapeHCurl<TSHAPE>::Shape(node,shapedata,phiHCurl,curlHCurl);

            // std::cout << "phiHCurl = " << phiHCurl << std::endl;
            // std::cout << "curlHCurl =  " << curlHCurl << std::endl;
            
            //Compute the curl for each edge
            top::ComputeConstantHCurl(node,N0Function,curlN0);

            // std::cout << "Constant phi = " << N0Function << std::endl;
            // std::cout << "Constant Curl = " << curlN0 << std::endl;
            // Checks if all edges have the same curl value.
            bool condHcurl = true;
            bool condHcurlN0 = true;
            for (int i = 0; i < dim; i++){
                for (int j = 0; j < nEdges; j++){
                    REAL funVal = (phiHCurl(i,2*j  )+phiHCurl(i,2*j+1)) / 2.;
                    if (fabs(funVal-N0Function(i,j)) > tol){
                        condHcurl = false;
                    } 
                }
            }
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < nEdges; j++){
                    REAL curlVal = (curlHCurl(i,2*j  )+curlHCurl(i,2*j+1))/2.;
                    if (fabs(curlVal-curlN0(i,j)) > tol){
                        condHcurlN0 = false;
                    } 
                }
            }
            if(!condHcurl || !condHcurlN0){
                std::cerr << "\n" + testName + " failed Hcurl" +
                                "\ntopology: " + MElementType_Name(type);
            }

            REQUIRE(condHcurl);
            REQUIRE(condHcurlN0);
        }
    }//Testing Constant Curl

}//namespace
