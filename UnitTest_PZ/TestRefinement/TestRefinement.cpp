//
//  TestRefinement.cpp
//  PZ
//
//  Created by Pedro Lima on 15/apr/2021.
//

// #include <pz_config.h>
#include "TPZRefPatternDataBase.h"
#include "pzgmesh.h"

#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

#include <stdio.h>
#include <iostream>
#define Zero -1.e-12;

// #ifdef PZ_USING_BOOST

#ifndef WIN32
    #define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MAIN pz refinement tests

#include <boost/test/unit_test.hpp>

using namespace pztopology;
// Forward declaring
namespace refinementtests{
    template <class topology>
    void TestingUniformRefinements(bool RefFromDatabase = 1);
}

BOOST_AUTO_TEST_SUITE(refinement_tests)

    BOOST_AUTO_TEST_CASE(refinement_tests_1)
    {
        // Test hardcoded refinements
        refinementtests::TestingUniformRefinements<TPZLine>(0);
        refinementtests::TestingUniformRefinements<TPZTriangle>(0);
        refinementtests::TestingUniformRefinements<TPZQuadrilateral>(0);
        refinementtests::TestingUniformRefinements<TPZTetrahedron>(0);
        refinementtests::TestingUniformRefinements<TPZCube>(0);
        refinementtests::TestingUniformRefinements<TPZPrism>(0);
        refinementtests::TestingUniformRefinements<TPZPyramid>(0);
        // Test refinements from database
        refinementtests::TestingUniformRefinements<TPZLine>(1);
        refinementtests::TestingUniformRefinements<TPZTriangle>(1);
        refinementtests::TestingUniformRefinements<TPZQuadrilateral>(1);
        refinementtests::TestingUniformRefinements<TPZTetrahedron>(1);
        refinementtests::TestingUniformRefinements<TPZCube>(1);
        refinementtests::TestingUniformRefinements<TPZPrism>(1);
        refinementtests::TestingUniformRefinements<TPZPyramid>(1);
    }


BOOST_AUTO_TEST_SUITE_END()

namespace refinementtests{

    void TestJacobian(TPZGeoMesh* gmesh);

    template<class TTopology>
    void TestingUniformRefinements(bool RefFromDatabase){

        TPZGeoMesh* gmesh = new TPZGeoMesh;
        // using TTopology pztopology::TTopology;

        // Create nodes
        int nnodes = TTopology::NCornerNodes;
        TPZManVector<REAL,3> coord(3,0.);
        TPZManVector<int64_t,6> cornerindices(nnodes,-1);
        for(int j=0; j < nnodes; j++){
            cornerindices[j] = j;
            gmesh->NodeVec().AllocateNewElement();
            TTopology::ParametricDomainNodeCoord(j,coord);
            coord.resize(3);
            gmesh->NodeVec()[j].Initialize(coord,*gmesh);
        }

        // Create father element
        gRefDBase.InitializeUniformRefPattern(TTopology::Type());
        int64_t index = 0;
        TPZGeoEl* gel = gmesh->CreateGeoElement(TTopology::Type(),cornerindices,1,index,RefFromDatabase);
        gmesh->BuildConnectivity();
        
        // Refine
        TPZManVector<TPZGeoEl*,10> children;
        gel->Divide(children);

        // Test subelements:

        // Check determinant of jacobian matrix
        TestJacobian(gmesh);
        
        /** Suggestions for future tests to add:
         * - Initialize and test ALL refinement patterns from the database
        */


        delete gmesh;
    }







    void TestJacobian(TPZGeoMesh* gmesh){

        const REAL notcomputed = -1024.0;
        static std::string testName = __PRETTY_FUNCTION__;

        for(auto gel : gmesh->ElementVec()){
            if(!gel) continue;
            if(gel->HasSubElement()) continue;

            TPZVec<REAL> qsi(gel->Dimension(),0.); TPZFMatrix<double> jac; TPZFMatrix<double> axes; REAL detjac = notcomputed; TPZFMatrix<double> jacinv;
            gel->Jacobian(qsi,jac,axes,detjac,jacinv);
            bool cond = detjac > Zero;
            BOOST_CHECK_MESSAGE(cond,
                                "\n"+testName+" failed"+
                                "\nUniform Refinement: "+gel->EldestAncestor()->TypeName()+"\n"+
                                "Subelement: "+ std::to_string(gel->WhichSubel())+" has negative detjacobian.\n"
                                );
        }
    }
}//namespace
// #endif
