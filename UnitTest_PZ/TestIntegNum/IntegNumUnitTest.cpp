/**
 * @file
 * @brief Contains Unit Test in Boost framework for Numerical integration module of the NeoPZ
 */

#include "pzvec.h"
#include "pzquad.h"

#include <iostream>
#include <fstream>
using namespace std;
#include <string>

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz numericintegration tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
// To compare test results with the stored values in output file with the expected results
#include "boost/test/output_test_stream.hpp"

#endif

#define MAXORDER 13

// Output file with the first test.
int FirstOrder = 0;
int max_order = MAXORDER;
int expo = 2;
int NDigitsPrec = 14;

std::string dirname = PZSOURCEDIR;

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
double Funcao(TPZVec<REAL> &point, int expo, int p) {
	REAL val = 0.0, power;
	REAL x = point[0];
	REAL y = point[1];
	REAL z = point[2];
	int ii, jj;
	// The degree p can not to be lower than expo
	if(p<expo) expo = p;
	val = CoefficientX(0,p);
	if(!p) return val;
	val += ((CoefficientX(1,p)*x)+(CoefficientY(1,p)*y)+(CoefficientZ(1,p)*z));
	if(p==1) {
		return val;
	}
	for(ii=2; ii<= expo; ii++) {
		val += (CoefficientX(ii,p) * x * x);
		val += (CoefficientY(ii,p) * y * y);
		val += (CoefficientZ(ii,p) * z * z);
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
void TestingNumericIntegrationRule(int exp, int p,boost::test_tools::output_test_stream &out) {
	
	TPZVec<REAL> point(3,0.);
	REAL weight = 0.;
	REAL integral = 0.;
	char text[64];
	char format[32];
	TPZVec<int> order(3,p);
	memset(text,0,strlen(text));
	memset(format,0,strlen(format));

	// Integration rule
	NumInteg intrule(p);
	intrule.SetOrder(order);

	int npoints = intrule.NPoints();
	int it;

	cout << "\nIntegration dimension = " << intrule.Dimension() << "\t order = " << p << "\t NPoints = " << npoints << "\n";
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
//		cout << "\t Point " << it << " : " << point[0] << "\t" << point[1] << "\t" << point[2] << std::endl;
//		cout << "\t Weight : " << weight << std::endl;
		integral += weight * Funcao(point,exp,p);
	}
	
	sprintf(format,"%%.%dlf",NDigitsPrec);
	sprintf(text,format,integral);
	cout << "Numerical integration value = " << text << std::endl;
	// Stores the obtained integral value
	out << text << "\n";

	BOOST_CHECK( out.match_pattern() );
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

	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult1D.txt";
	boost::test_tools::output_test_stream output(filename,true);
	
	int order;
	// Whether FirstResult.txt is empty then it is a first run
	for(order=FirstOrder;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZInt1d>(expo,order,output);   // OK
	}
	// Conclusion: It's failed at order = 1. But the erro is proportional e-10 from order = 8
}

BOOST_AUTO_TEST_CASE(numinteg2D_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult2D.txt";
	boost::test_tools::output_test_stream output(filename,true);
	
	int order;
	// Quadrilateral parametric space
	for(order=FirstOrder;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZIntQuad>(expo,order,output);   // OK
	}
	// Conclusion: We have problem at order 1. But the erro is proportional e-10 from order = 8
	output << "\n";
	BOOST_CHECK( output.match_pattern() );	
	
	// Triangular parametric space
	for(order=FirstOrder;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZIntTriang>(expo,order,output);
	}
	// Conclusion: We have problem at order 1 and 5. It's failed at order >= 10.
}

BOOST_AUTO_TEST_CASE(numinteg3D_tests) {

	if(NDigitsPrec > 9)
		NDigitsPrec = 9;
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3D.txt";
	boost::test_tools::output_test_stream output(filename,true);
	
	int order;
	// Cube
	for(order=FirstOrder;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntCube3D>(expo,order,output);   // OK
	}
	// Conclusion: We have problem at order 1 and 0. But the erro is proportional e-10 from order = 8
	output << "\n";
	BOOST_CHECK( output.match_pattern() );
	
	// Tetrahedram
	for(order=FirstOrder;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntTetra3D>(expo,order,output);
	}
	// Conclusion: It's failed at order = 1 and order = 4
	output << "\n";	
	BOOST_CHECK( output.match_pattern() );
	
	// Pyramidal
	for(order=FirstOrder;order<max_order;order++) {
		//		TestingNumericIntegrationRule<TPZIntPyram3D>(expo,order,out);  // It isn't implemented yet
	}
	
	// Prism
	for(order=FirstOrder;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntPrism3D>(expo,order,output);
	}
	// Conclusion: It's failed at order = 1 and order = 0
}

BOOST_AUTO_TEST_SUITE_END()

#endif
