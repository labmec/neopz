/**
 * @file
 * @brief Unit test for numerical integration module.
 */

#include "pzvec.h"
#include "pzquad.h"

//IO
#include <iostream>
using namespace std;

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz numericintegration tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#endif

#define MAXORDER 13

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
	NumInteg intrule(order);
	int npoints = intrule.NPoints();
	int it;
	cout << "\nIntegration dimension = " << intrule.Dimension() << "\torder = " << order << "\tNPoints = " << npoints << "\n";
	for (it=0;it<npoints;it++){
		intrule.Point(it,point,weight);
		cout << "\tPoint " << it << " : " << point[0] << "\t" << point[1] << "\t" << point[2] << std::endl;
		cout << "\t Weight : " << weight << std::endl;
		integral += weight * Funcao(point,exp,order);
	}
	
	cout << "Numerical integration value = " << integral << "\tDiff = " << fabs(intvalue - integral) << std::endl << std::endl;

	BOOST_CHECK(IsZero(intvalue-integral));
}

/**
 * Using Mathematica, I had the following values:
 * (1D) i=2, p=6, from x = -1 until x = 1 then  integral = 2.169215575174690848712687509073837669906381782925
 * para 2D:
 * i=2, p=6, from x = -1 until x = 1, and y=-1 until y=1 then  integral = 3.6230005930154684018715427138244729500584455281477 (quadrilateral)
 * i=2, p=6, from x = 0. until x = 1, and y=0 until y=1-x then integral = 0.2645718117897931699940052644354916134028121112478236818865112770629 (triangular)
 * para 3D:
 * i=2, p=6, from x = -1 until x = 1, and y=-1 until y=1 and z=-1 until z=1 then  integral = 9.6707943242327546581324717367469321473229043134898 (cube)
 */

BOOST_AUTO_TEST_SUITE(numinteg_tests)

BOOST_AUTO_TEST_CASE(numinteg1D_tests) {
	int i = 2;
	int order;
	REAL intvalue[MAXORDER];
	intvalue[0] = 2.000000000000000;
	intvalue[1] = 2.606198284550454;
	intvalue[2] = 2.659572164415588;
	intvalue[3] = 1.499341835485549;
	intvalue[4] = 2.367617787494460;
	intvalue[5] = 1.825083430864047;
	intvalue[6] = 2.169215575174691;
	intvalue[7] = 1.617745418673051;
	intvalue[8] = 2.480691807001154;
	intvalue[9] = 1.347699766137747;
	intvalue[10] = 1.417801801857337;
	intvalue[11] = 1.935192061654517;
	intvalue[12] = 1.429663752832786;
	for(order=0;order<MAXORDER;order++) {
		TestingNumericIntegrationRule<TPZInt1d>(i,order,intvalue[order]);   // OK
	}
	// Conclusion: It's failed at order = 1. But the erro is proportional e-10 from order = 8
}

BOOST_AUTO_TEST_CASE(numinteg2D_tests) {
	int i = 2;
	int order;
	// Quadrilateral parametric space
	REAL intvalue[MAXORDER];
	intvalue[0] = 4.;
	intvalue[1] = 6.424793138201817;
	intvalue[2] = 4.310074335087271;
	intvalue[3] = 2.626129673372531;
	intvalue[4] = 6.054379903820096;
	intvalue[5] = 2.924805380542269;
	intvalue[6] = 3.6230005930154684018715427138244729500584455281477;
	intvalue[7] = 4.556300644939264;
	intvalue[8] = 4.577512525115554;
	intvalue[9] = 1.6940832032465925;
	intvalue[10] = 4.052863938018176;
	intvalue[11] = 3.858582377588495;
	intvalue[12] = 1.6518896896567403;

	for(order=0;order<MAXORDER;order++) {
		TestingNumericIntegrationRule<TPZIntQuad>(i,order,intvalue[order]);   // OK
	}
	// Conclusion: We have problem at order 1. But the erro is proportional e-10 from order = 8

	// Triangular parametric space
	intvalue[0] = 0.5;
	intvalue[1] = 0.9320398994069123;
	intvalue[2] = 0.54479546786258;
	intvalue[3] = 0.5063395201360538;
	intvalue[4] = 0.45428110865992366;
	intvalue[5] = 0.2509209988237399;
	intvalue[6] = 0.26457181178979317;
	intvalue[7] = 0.4853077813019233;
	intvalue[8] = 0.8543252469564175;
	intvalue[9] = 0.3195852820307958;
	intvalue[10] = 0.3282395374596146;
	intvalue[11] = 0.49096556829476024;
	intvalue[12] = 0.18197735362039558;
	for(order=0;order<MAXORDER;order++) {
		TestingNumericIntegrationRule<TPZIntTriang>(i,order,intvalue[order]);
	}
	// Conclusion: We have problem at order 1 and 5. It's failed at order >= 10.
}

BOOST_AUTO_TEST_CASE(numinteg3D_tests) {
	int i = 2;
	int order;
	// Cube
	REAL intvalue[MAXORDER];
	intvalue[0] = 10.42479313820182;
	intvalue[1] = 15.27437941460545;
	intvalue[2] = 11.04494180837636;
	intvalue[3] = 7.677052484946879;
	intvalue[4] = 14.53355294584201;
	intvalue[5] = 8.274403899286355;
	intvalue[6] = 9.670794324232754;
	intvalue[7] = 11.53739442808034;
	intvalue[8] = 11.57981818843293;
	intvalue[9] = 5.812959544695003;
	intvalue[10] = 10.53052101423817;
	intvalue[11] = 10.14195789337881;
	intvalue[12] = 5.728572517515299;
	for(order=0;order<MAXORDER;order++) {
		TestingNumericIntegrationRule<TPZIntCube3D>(i,order,intvalue[order]);   // OK
	}
	// Conclusion: We have problem at order 1. But the erro is proportional e-10 from order = 8

	// Tetrahedram
	intvalue[0] = 0.21688291481409032598;
	intvalue[1] = 0.31731541110893806824;
	intvalue[2] = 0.2271127993891759;
	intvalue[3] = 0.2227611387537444;
	intvalue[4] = 0.2227611387537444;
	intvalue[5] = 0.1579731476556451;
	intvalue[6] = 0.1592039908614216;
	intvalue[7] = 0.2114714205797953;
	intvalue[8] = 0.3036594998229580;
	intvalue[9] = 0.1789852244559073;
	intvalue[10] = 0.1737775985061321;
	intvalue[11] = 0.2150662368602839;
	intvalue[12] = 0.1447150978897556;
	for(order=0;order<MAXORDER;order++)
		TestingNumericIntegrationRule<TPZIntTetra3D>(i,order,intvalue[order]);
	// Conclusion: It's failed at order = 1 and order = 4

/*
	// Pyramidal
	intvalue[0] = ;
	intvalue[1] = ;
	intvalue[2] = ;
	intvalue[3] = ;
	intvalue[4] = ;
	intvalue[5] = ;
	intvalue[6] = ;
	intvalue[7] = ;
	intvalue[8] = ;
	intvalue[9] = ;
	intvalue[10] = ;
	intvalue[11] = ;
	intvalue[12] = ;
	for(order=0;order<MAXORDER;order++)
		TestingNumericIntegrationRule<TPZIntPyram3D>(i,order,intvalue[order]);
*/	
	// Prism
	intvalue[0] = 1.303099142275227;
	intvalue[1] = 2.167178941089052;
	intvalue[2] = 1.392690078000387;
	intvalue[3] = 1.315778182547335;
	intvalue[4] = 1.211661359595075;
	intvalue[5] = 0.8049411399227070;
	intvalue[6] = 0.8322427658548136;
	intvalue[7] = 1.273714704879074;
	intvalue[8] = 2.011749636188062;
	intvalue[9] = 0.9422697063368188;
	intvalue[10] = 0.9595782171944565;
	intvalue[11] = 1.285030278864748;
	intvalue[12] = 0.6670538495160184;
	for(order=0;order<MAXORDER;order++)
		TestingNumericIntegrationRule<TPZIntPrism3D>(i,order,intvalue[order]);
	// Conclusion: It's failed at order = 1 and order = 0
}

BOOST_AUTO_TEST_SUITE_END()
	
#endif