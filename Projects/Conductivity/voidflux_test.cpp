//
//  voidflux_test.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz voidflux tests

//#include <boost/test/unit_test.hpp>

//#include <boost/test/output/compiler_log_formatter.hpp>

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "xcode_log_formatter.hpp"

using namespace boost::unit_test;

struct xcode_config
{
	xcode_config()
	{
		unit_test_log.set_formatter( new xcode_log_formatter );
	}
	
	~xcode_config() {}
};

BOOST_GLOBAL_FIXTURE(xcode_config);


#include "TPZConductivityProblem.h"
#include "pzgmesh.h"
#include <iostream>


// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(conductivity)

BOOST_AUTO_TEST_CASE(geometry_test)
{
    TPZConductivityProblem problem;
    problem.SetDomainSize(1.,1.);
    problem.SetMesh(5,5);
    TPZAutoPointer<TPZCompMesh> cmesh = problem.GenerateCompMesh();
    TPZGeoMesh * gmesh = cmesh->Reference();
    REAL perimeter = problem.Perimeter(*gmesh);
    REAL area = problem.DomainArea(*gmesh);
    BOOST_CHECK_SMALL(area-(REAL)(1.L), (REAL)1.e-6L);
    BOOST_CHECK_SMALL(perimeter-(REAL)(20.L), (REAL)1.e-6L);
}

BOOST_AUTO_TEST_CASE(linearity_test)
{
    TPZConductivityProblem problem;
    TPZManVector<REAL,3> delx;
    problem.GetDomainSize(delx);
    TPZManVector<int,3> nx;
    problem.GetMesh(nx);
    REAL flux = problem.ComputeFlux();
    delx[1] *= 2.;
    problem.SetDomainSize(delx[0],delx[1]);
    REAL flux2 = problem.ComputeFlux();
    BOOST_CHECK_SMALL((flux- flux2/((REAL)(2.L)))/flux, (REAL)1.e-6L);
    REAL conductivity = problem.GetConductivity();
    conductivity *= 2.;
    problem.SetConductivity(conductivity);
    REAL flux3 = problem.ComputeFlux();
    BOOST_CHECK_SMALL((flux2- flux3/((REAL)(2.0L)))/flux2, (REAL)1.e-6L);
}


BOOST_AUTO_TEST_SUITE_END()


#else

#include <iostream>

int main()
{
    std::cout << "Boost needs to be configured for testing to work\n";
    return 0;
}

#endif