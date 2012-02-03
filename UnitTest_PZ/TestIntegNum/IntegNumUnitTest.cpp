/**
 * @file
 * @brief Contains Unit Test in Boost framework for Numerical integration module of the NeoPZ
 */

#include "pzvec.h"
#include "pzquad.h"

#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
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

// Global variables to test
int max_order = MAXORDER;
int expo = 2;
unsigned int NDigitsPrec = 18;

std::string dirname = PZSOURCEDIR;


// Computing coefficient of the polinomial function a(i,p) = sin(i*p*p)
REAL CoefficientX(int i, int p) {
	if(!i) return (REAL)(0.125L);
	return (REAL)(i*p*p);
}
REAL CoefficientY(int i, int p) {
	return (REAL)(i*p);
}
REAL CoefficientZ(int i, int p) {
	return (REAL)(i);
}

// Defining a function to integrate f(x,y,z,i,p) = CoeffX(0,p) + CoeffX(1,p)x + CoeffX(2,p)x^2 + ... + CoeffX(i,p)x^i 
//                                                             + CoeffY(1,p)y + CoeffY(2,p)y^2 + ... + CoeffY(i,p)y^i
//                                                             + CoeffZ(1,p)z + CoeffZ(2,p)z^2 + ... + CoeffZ(i,p)z^i 
REAL Funcao(TPZVec<REAL> &point, int expo, int p) {
	REAL val = (REAL)(0.0L);
	REAL x = point[0];
	REAL y = point[1];
	REAL z = point[2];
	int ii;
	// The degree p can not to be lower than expo
	if(p<expo) expo = p;
	val = CoefficientX(0,p);
	if(!p) return val;
	val = val + ((CoefficientX(1,p)*x)+(CoefficientY(1,p)*y)+(CoefficientZ(1,p)*z));
	if(p==1) {
		return val;
	}
	for(ii=2; ii<= expo; ii++) {
		val = val + (CoefficientX(ii,p) * x * x);
		val = val + (CoefficientY(ii,p) * y * y);
		val = val + (CoefficientZ(ii,p) * z * z);
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
	// Variables to computing numerical integration
	TPZVec<REAL> point(3,0.L);
	REAL weight = 0.L;
	REAL integral = 0.L;
	TPZVec<int> order(3,p);

	// Variables to put numerical value as wish format
	char text[64];
	char format[32];
	memset(text,0,strlen(text));
	memset(format,0,strlen(format));
	
	// Integration rule
	NumInteg intrule(p);
	intrule.SetOrder(order);

	int npoints = intrule.NPoints();
	int it;
	
	// Computing the definite integral for Funcao on parametric element
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
		integral = integral + weight * Funcao(point,exp,p);
	}

	// USING TEXT FORMATED TO CHECK WITH FILE WITH PREVIOUS INTEGRATION DATA
	// Formating the integral value obtained
	if(integral < 10.) sprintf(format,"%%.%dLf",NDigitsPrec);
	else if(integral < 100.) sprintf(format,"%%.%dLf",NDigitsPrec-1);
	else if(integral < 1000.) sprintf(format,"%%.%dLf",NDigitsPrec-2);
	else if(integral < 10000.) sprintf(format,"%%.%dLf",NDigitsPrec-3);
	else sprintf(format,"%%.%dLf",NDigitsPrec-4);
	sprintf(text,format,integral);
	
	// Stores the obtained integral value into the boost::output to compare with match_pattern()
	out << text << "\n";

	BOOST_CHECK_MESSAGE( out.match_pattern() , "\nIntegration: Dim = " << intrule.Dimension() << "\t Order = " << p << "\t NPoints = " << npoints << "\t Value = " << text << "\n");
}
template <class NumInteg>
void TestingNumericIntegrationRule(int exp, int p,std::ifstream &input) {
	// Variables to computing numerical integration
	TPZVec<REAL> point(3,0.L);
	REAL weight = 0.L;
	REAL integral = 0.L;
	TPZVec<int> order(3,p);
		
	// Integration rule
	NumInteg intrule(p);
	intrule.SetOrder(order);
	
	int npoints = intrule.NPoints();
	int it;
	
	// Computing the definite integral for Funcao on parametric element
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
		integral = integral + weight * Funcao(point,exp,p);
	}

	// USING IMPORTED DATA FROM FILE WITH WISHED VALUES
	// Variables to import data from file with wished data
	REAL inputdata;
	long double tol = 1.e-17L;
	input >> inputdata;
	// Making tol compatible with the wished significant digits
	for(it=0;it<(18-NDigitsPrec);it++)
		tol *= 10.L;
	if(inputdata > 10.0)
		tol *= 10.L;
	if(inputdata > 100.0)
		tol *= 10.L;
	if(inputdata > 1000.0)
		tol *= 10.L;

	// SEE -> check the predicate, if the predicate is false then the message will be displayed.
	BOOST_CHECK_MESSAGE( fabs(integral-inputdata) < tol , "\nIntegration: Dim = " << intrule.Dimension() << "\t Order = " << p << "\t NPoints = " << npoints << "\t Value = " << integral << " - " << fabs(integral-inputdata) << "\n");
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

	// File with integration values calculated previously
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult1D.txt";
	std::ifstream olddata(filename.c_str());

	// Whether FirstResult.txt is empty then it is a first run
	for(int order=0;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZInt1d>(expo,order,olddata);   // OK
	}
}

BOOST_AUTO_TEST_CASE(numinteg2D_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult2D.txt";
	std::ifstream olddata(filename.c_str());

	int order;
	// Quadrilateral parametric space
	for(order=0;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZIntQuad>(expo,order,olddata);   // OK
	}
	
	// Triangular parametric space
	NDigitsPrec -= 7;                       // because the precision of the points and weights is to 15 digits
	for(order=0;order<max_order-2;order++) {
		TestingNumericIntegrationRule<TPZIntTriang>(expo,order,olddata);   // OK
	}
	NDigitsPrec += 7;
}

BOOST_AUTO_TEST_CASE(numinteg3D_tests) {

	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3D.txt";
	std::ifstream olddata(filename.c_str());

	int order;
	// Cube
	for(order=0;order<max_order;order++) {
		if(!order) {
			continue;
		}
		if(order==12) NDigitsPrec -= 2;
		TestingNumericIntegrationRule<TPZIntCube3D>(expo,order,olddata);   // OK
		if(order==12) NDigitsPrec += 2;
	}
	
	// Tetrahedram
	NDigitsPrec -= 9;                       // because the precision of the points and weights is to 15 digits
	for(order=0;order<max_order;order++) {
		if(!order) {
			continue;
		}
		TestingNumericIntegrationRule<TPZIntTetra3D>(expo,order,olddata);   // OK
	}
	NDigitsPrec += 9;
	
	// Pyramidal
	NDigitsPrec -= 9;                       // because the precision of the points and weights is to 15 digits
	for(order=2;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZIntPyram3D>(expo,order,olddata);   // OK
	}
	NDigitsPrec += 9;
	
	// Prism
	NDigitsPrec -= 8;                       // because the precision of the points and weights is to 15 digits
	for(order=0;order<max_order-2;order++) {
		if(!order) {
			continue;
		}
		TestingNumericIntegrationRule<TPZIntPrism3D>(expo,order,olddata);   // OK
	}
	NDigitsPrec += 8;
}

BOOST_AUTO_TEST_SUITE_END()

/**
 
BOOST_AUTO_TEST_CASE(numinteg1DM_tests) {
	
	// 
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult1D.txt";
	boost::test_tools::output_test_stream output(filename,true);
	
	// Whether FirstResult.txt is empty then it is a first run
	for(int order=0;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZInt1d>(expo,order,output);
	}
}

BOOST_AUTO_TEST_CASE(numinteg2DM_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult2D.txt";
	boost::test_tools::output_test_stream output(filename,true);
	
	int order;
	// Quadrilateral parametric space
	for(order=0;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZIntQuad>(expo,order,output);
	}
	output << "\n";
	BOOST_CHECK( output.match_pattern() );	
	
	// Triangular parametric space
	for(order=0;order<max_order;order++) {
		TestingNumericIntegrationRule<TPZIntTriang>(expo,order,output);
	}
}

BOOST_AUTO_TEST_CASE(numinteg3DM_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3D.txt";
	boost::test_tools::output_test_stream output(filename,true);
	
	int order;
	// Cube
	for(order=0;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntCube3D>(expo,order,output);
	}
	output << "\n";
	BOOST_CHECK( output.match_pattern() );
	
	// Tetrahedram
	for(order=0;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntTetra3D>(expo,order,output);
	}
	output << "\n";	
	BOOST_CHECK( output.match_pattern() );
	
	// Pyramidal
	for(order=0;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntPyram3D>(expo,order,output);
	}
	output << "\n";	
	BOOST_CHECK( output.match_pattern() );
	
	// Prism
	for(order=0;order<max_order;order++) {
		if(!order) {
			output << "\n";
			BOOST_CHECK( output.match_pattern() );	
			continue;
		}
		TestingNumericIntegrationRule<TPZIntPrism3D>(expo,order,output);
	}
}
*/

#endif


// Compute the points and weights for Gauss Legendre Quadrature over the parametric 1D element [-1.0, 1.0]
void ComputingGaussLegendreQuadrature(int Npoints,ofstream &GaussLegQuadrature) {
	const long double tol = 1.e-16;
	long double z1, z, pp, p3, p2, p1;
	int i;
	TPZVec<long double> Points(Npoints);
	TPZVec<long double> Weights(Npoints);
	
	int m = (Npoints+1)*0.5;
	
	for(i=0;i<m;i++) {
		p1 = ((REAL)i)+(REAL)0.75;
		p2 = ((REAL)Npoints)+(REAL)(0.5);
		
		z = cosl(((REAL)M_PI)*p1/p2);
		do {
			p1 = (REAL)1.0;
			p2 = (REAL)0.0;
			for(int j=0;j<Npoints;j++) {
				p3 = p2;
				p2 = p1;
				p1 = (((REAL)(2.0*j+1.0))*z*p2-((REAL)j)*p3)/((REAL)(j+1));
			}
			pp = ((REAL)Npoints)*(z*p1-p2)/(z*z-((REAL)1.0));
			z1 = z;
			z = z1-p1/pp;
		}while(fabs(z-z1) > tol);
		Points[i] = -z;
		Points[Npoints-1-i] = z;
		Weights[i] = 2.0/((1.0-z*z)*pp*pp);
		Weights[Npoints-1-i] = Weights[i];
	}
	// Printing points and weights
	char text[64];
	char format[32];
	memset(text,0,strlen(text));
	memset(format,0,strlen(format));
	
	GaussLegQuadrature << Npoints << "\n";
	sprintf(format,"%%.%dLf",25);
	for(i=0;i<Npoints;i++) {
		sprintf(text,format,Points[i]);
		GaussLegQuadrature << text << "\t";
		sprintf(text,format,Weights[i]);
		GaussLegQuadrature << text << "\n";
	}
	GaussLegQuadrature << "\n";
}
