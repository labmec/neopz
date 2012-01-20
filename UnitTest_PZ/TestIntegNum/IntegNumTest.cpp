/*
 *  IntegNumTest.cpp
 *  PZ
 *
 *  Created by Jorge on 1/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

//#include "IntegNumTest.h"
#include "pzvec.h"
#include "pzquad.h"

//IO
#include <iostream>
using namespace std;

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz numericintegration tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

// Computing coefficient of the polinomial function a(i,p) = sin(i*p*p)
double CoefficientX(int i, int p) {
	if(!i) return 1.0;
	return sin((double)(i*p*p));
}
double CoefficientY(int i, int p) {
	return sin((double)(i*p));
}
double CoefficientZ(int i, int p) {
	return sin((double)(i));
}
// Defining a function to integrate f(x,y,z,i,p) = CoeffX(0,p) + CoeffX(1,p)x + CoeffX(2,p)x^2 + ... + CoeffX(i,p)x^i 
//                                                             + CoeffY(1,p)y + CoeffY(2,p)y^2 + ... + CoeffY(i,p)y^i
//                                                             + CoeffZ(1,p)z + CoeffZ(2,p)z^2 + ... + CoeffZ(i,p)z^i 
double Funcao(TPZVec<REAL> &point, int i, int p) {
	REAL val = 0.0;
	REAL x = point[0];
	REAL y = point[1];
	REAL z = point[2];
	val = CoefficientX(0,p);
	for(int ii=1; ii<= i; ii++) {
		val += (CoefficientX(ii,p) * pow(x,ii));
		val += (CoefficientY(ii,p) * pow(y,ii));
		val += (CoefficientZ(ii,p) * pow(z,ii));
	}
	return val;
}

#ifdef USING_BOOST

/**
 * @brief Tests the numerical integration rule implemented in neopz.lib. 
 * @param order Order of the numeric integration
 * @note The argument of the template is the Numeric Integration Rule.
 */
/**
 * Integrando numericamente a funcao f(x,y,z) = a(0,p) + a(1,p)x + a(2,p)x*x + ... + a(i,p)x^i + a(1,p)y + a(2,p)y*y + ...
 *                                        + a(i,p)y^i + a(1,p)z + a(2,p)z^2 + ... + a(i,p)z^i 
 * O matemathica fornece os seguintes valores:
 * Para i = 2, order = 6, 1D, integral = 
 * Para i = 2, order = 6, 2D, integral = 
 */
template <class NumInteg>
void TestingNumericIntegrationRule(int exp, int order, REAL intvalue) {
	TPZVec<REAL> point(3,0.);
	REAL weight = 0.;
	REAL integral = 0.;
	
	//=====1D Rule=====================================
	NumInteg ordem1d (order);
	int npoints = ordem1d.NPoints();
	int it;
	for (it=0;it<npoints;it++){
		ordem1d.Point(it,point,weight);
		integral += weight * Funcao(point,exp,order);
	}
	
	cout << "Numerical integration value = " << integral << endl;
	cout << "Integration order = " << order << "\tNPoints = " << npoints << "\n";
	for(it=0;it<npoints;it++) {
		cout << "\tPoint " << it << " : " << point[0] << "\t" << point[1] << "\t" << point[2] << endl;
		cout << "\t Weight : " << weight << endl;
	}
	cout << "\t Diferença = " << fabs(integral - intvalue) << endl; 

	BOOST_CHECK(IsZero(intvalue-integral));
}

/**
 * Using Mathematica, we had the following values:
 * i=2, p=6, from x = -1 until x = 1 then  integral = 2.169215575174690848712687509073837669906381782925
 * para 2D:
 * i=2, p=6, from x = -1 until x = 1, and y=-1 until y=1 then  integral = 3.6230005930154684018715427138244729500584455281477
 * para 3D:
 * i=2, p=6, from x = -1 until x = 1, and y=-1 until y=1 and z=-1 until z=1 then  integral = 9.6707943242327546581324717367469321473229043134898 
 */

BOOST_AUTO_TEST_SUITE(numinteg_tests)

BOOST_AUTO_TEST_CASE(numinteg1D_tests) {
	int i = 2;
	int order = 6;
	REAL intvalue = 2.169215575174690848712687509073837669906381782925;
	TestingNumericIntegrationRule<TPZInt1d>(i,order,intvalue);   // OK
}

BOOST_AUTO_TEST_CASE(numinteg2D_tests) {
	int i = 2;
	int order = 6;
	REAL intvalue = 3.6230005930154684018715427138244729500584455281477;
	TestingNumericIntegrationRule<TPZIntTriang>(i,order,intvalue);
	TestingNumericIntegrationRule<TPZIntQuad>(i,order,intvalue);   // OK
}

BOOST_AUTO_TEST_CASE(numinteg3D_tests) {
	int i = 2;
	int order = 6;
	REAL intvalue = 9.6707943242327546581324717367469321473229043134898;
	TestingNumericIntegrationRule<TPZIntTetra3D>(i,order,intvalue);
	TestingNumericIntegrationRule<TPZIntPyram3D>(i,order,intvalue);
	TestingNumericIntegrationRule<TPZIntPrism3D>(i,order,intvalue);
	TestingNumericIntegrationRule<TPZIntCube3D>(i,order,intvalue);   // OK
}


BOOST_AUTO_TEST_SUITE_END()
	
#endif