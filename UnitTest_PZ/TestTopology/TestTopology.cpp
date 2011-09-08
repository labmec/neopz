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

#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz topology tests

#include <boost/test/unit_test.hpp>

using namespace pztopology;

template <class top>
void TestingSideProjections2() {
    int dim = top::Dimension;
    int nsides = top::NSides;
    TPZManVector<REAL,3> pt1,pt2,pt3;
    int is;
    for (is=0; is<nsides; is++) {
        int sidedim = top::SideDimension(is);
        TPZManVector<REAL,3> pt1(dim),pt2(dim),pt3(sidedim);
        top::CenterPoint(nsides-1,pt1);
        TPZTransform tr1 = top::TransformSideToElement(is);
        TPZTransform tr2 = top::TransformElementToSide(is);
        TPZTransform tr3;
        tr3 = tr2.Multiply(tr1);
        BOOST_CHECK_EQUAL(tr3.Mult().Rows(), sidedim);
    }
}

// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(topology_tests2)

BOOST_AUTO_TEST_CASE(projection_tests2)
{
    TestingSideProjections2<TPZLine>();
    // Ensure that subtracting 6 from 8 gives 2.
    //TPZVoidFlux flux(1,2.);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
