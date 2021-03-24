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

#ifdef PZ_USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz topology tests

#include <boost/test/unit_test.hpp>

using namespace pztopology;
namespace topologytests{
    template <class top>
    void TestingSideProjections();
    template <class top>
    void TestingSideNodeProjections();
}

// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(topology_tests)

BOOST_AUTO_TEST_CASE(projection_tests_1)
{
    topologytests::TestingSideProjections<TPZLine>();
    topologytests::TestingSideProjections<TPZTriangle>();
    topologytests::TestingSideProjections<TPZQuadrilateral>();
    topologytests::TestingSideProjections<TPZTetrahedron>();
    topologytests::TestingSideProjections<TPZCube>();
    topologytests::TestingSideProjections<TPZPrism>();
    topologytests::TestingSideProjections<TPZPyramid>();
}

    BOOST_AUTO_TEST_CASE(projection_tests_2)
    {
        topologytests::TestingSideNodeProjections<TPZTriangle>();
        topologytests::TestingSideNodeProjections<TPZQuadrilateral>();
        topologytests::TestingSideNodeProjections<TPZTetrahedron>();
        topologytests::TestingSideNodeProjections<TPZCube>();
        topologytests::TestingSideNodeProjections<TPZPrism>();
        topologytests::TestingSideNodeProjections<TPZPyramid>();
    }

BOOST_AUTO_TEST_SUITE_END()

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
                    BOOST_CHECK_MESSAGE(cond,
                                        "\n"+testName+" failed"+
                                        "\ntopology: "+MElementType_Name(type)+"\n"+
                                        "side: "+ std::to_string(iSide)
                    );
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
                BOOST_CHECK_MESSAGE(cond,"\n"+testName+" failed"+
                                         "\ntopology: "+MElementType_Name(type)+"\n"+
                                         "side: "+ std::to_string(iSide)
                                         );
            }
        }
    }//TestingSideNodeProjections
}//namespace
#endif
