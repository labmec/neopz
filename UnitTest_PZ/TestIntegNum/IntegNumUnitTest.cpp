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

// Global variables to test
int expo = 2;
unsigned int NDigitsPrec = 13;
// Conclusion:	Using NDigitsPrec = 4 can to run with REAL = float
//				Using NDigitsPrec = 13 can to run with REAL = double
//				Using NDigitsPrec = 13 can to run with REAL = long double
//  FOR ALL CUBATURE RULES

std::string dirname = PZSOURCEDIR;

int NTypes = 2;

/** 
 * @name Testing a numeric integration rule of order p for function with argument and coefficients
 * depending on p.
 * @{ 
 */

REAL CoefficientX(int i, int p);
REAL CoefficientY(int i, int p);
REAL CoefficientZ(int i, int p);

/**
 * @brief Defining a function to integrate \f$ f(x,y,z,i,p) = Cx(0,p) + Cx(1,p)x + Cx(2,p)x^2 + ... + Cx(i,p)x^i \f$
 *      \f$ + Cy(1,p)y + Cy(2,p)y^2 + ... + Cy(i,p)y^i + Cz(1,p)z + Cz(2,p)z^2 + ... + Cz(i,p)z^i \f$
 */
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

/** @brief Compute the coefficient of the polinomial function \f$Cx(i,p) = i*p*p\f$ but \f$ Cx(0,p) = 0.125 \f$*/
REAL CoefficientX(int i, int p) {
	if(!i) return (REAL)(0.125L);
	return (REAL)(i*p*p);
}
/** @brief Compute the coefficient of the polinomial function \f$Cy(i,p) = i*p\f$ */
REAL CoefficientY(int i, int p) {
	return (REAL)(i*p);
}
/** @brief Compute the coefficient of the polinomial function \f$Cz(i,p) = i \f$ */
REAL CoefficientZ(int i, int p) {
	return (REAL)(i);
}

/* @} */

/** 
 * @name Testing a numeric integration rule of order p for polinomial of order N
 * @{ 
 */

/** @brief Compute the value of polinomial \f$ P(x) = (N+1)*x^N + N^2 \f$. \n
 * The integral over \f$ [a,b] \f$ is \f$ \int[P(x)dx]\sub{a} = (b^{N+1} - a^{N+1}) + N^2 (b-a) \f$
 */
REAL PolinomialN(int N,REAL x) {
	REAL valuex = 1.0L;
	// Computes the value of \f$ x^N \f$
	for(int i=0;i<N;i++)
		valuex *= x;
	return (((REAL)(N+1))*valuex + ((REAL)(N*N)));
}
/**
 * @brief Compute the integral value of PolinomialN over the interval \f$[LimInf=a, LimSup=b]\f$ 
 * @param N Degree of the polinomial
 */
REAL IntegralOfPolinomialN(int N,REAL LimInf, REAL LimSup) {
	REAL A = 1.0L, B = 1.0L;
	for(int i=0;i<(N+1);i++) {
		A *= LimInf;
		B *= LimSup;
	}
	return ((B-A)+(((REAL)(N*N))*(LimSup-LimInf)));
}
/**
 * @brief Linear transformation to transform a interval [-1.,1.] to [a,b] 
 * @param point Point within the master element [-1.,1.]
 * @param LimInf Initial point of the interval image of the master
 * @param LimSup Final point of the interval image of the master
 */
REAL LinearTransformFromMaster(REAL point,REAL LimInf, REAL LimSup) {
	return (0.5L * ((LimSup - LimInf)*point + (LimSup + LimInf)));
}
/**
 * @brief Linear transformation to transform a interval [-1.,1.] to [a,b] 
 * @param point Point within the master element [-1.,1.]
 * @param LimInf Initial point of the interval image of the master
 * @param LimSup Final point of the interval image of the master
 */
REAL DifferenceOfLinearTransformFromMaster(REAL point,REAL LimInf, REAL LimSup) {
	return (0.5L * (LimSup - LimInf));
}

/** @} */

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
void TestingNumericIntegrationRule(int exp, int p,int type,boost::test_tools::output_test_stream &out) {
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
	intrule.SetType(type);
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
void TestingNumericIntegrationRule(int exp, int p,int type,std::ifstream &input) {
	// Variables to computing numerical integration
	TPZVec<REAL> point(3,0.L);
	REAL weight = 0.L;
	REAL integral = 0.L;
	TPZVec<int> order(3,p);
		
	// Integration rule
	NumInteg intrule(p);
	intrule.SetOrder(order,type);
	
	unsigned int npoints = intrule.NPoints();
	unsigned int it;
	std::string namerule;
	intrule.Name(namerule);
	
	// Computing the definite integral for Funcao on parametric element
	cout << " Cubature rule: " << namerule << "  Type " << type << " \tOrder " << p << " \tNPoints " << npoints << std::endl;
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
	//	cout << "\tPoint : " << point[0] << "\t" << point[1] << "\t" << point[2] << "\t Weight " << weight << std::endl;
		integral = integral + weight * Funcao(point,exp,p);
	}

	// USING IMPORTED DATA FROM FILE WITH WISHED VALUES
	// Variables to import data from file with wished data
	long double inputdata;
	long double tol = 1.L;
	input >> inputdata;
	// Making tol compatible with the wished significant digits
	for(it=0;it<NDigitsPrec;it++)
		tol *= 0.1L;
	
	if(inputdata > 10.0)
		tol *= 10.L;
	if(inputdata > 100.0)
		tol *= 10.L;
	if(inputdata > 1000.0)
		tol *= 10.L;

	// SEE -> check the predicate, if the predicate is false then the message will be displayed.
	BOOST_CHECK_MESSAGE( fabs(integral-inputdata) < tol , "\nIntegration: Dim = " << intrule.Dimension() << "\t Order = " << p << "\t NPoints = " << npoints << "\t Value = " << integral << " - (" << fabs(integral-inputdata) << ")\n");
}

template<class NumInteg>
void ComparePointsWithNewRule(int p) {
	// Variables to computing numerical integration
	TPZVec<REAL> point(3,0.L);
	REAL weight = 0.L;
	TPZVec<int> order(3,p);
	long double tol = 1.0e-18;
	
	// Integration rule
	NumInteg intrule(p);
	intrule.SetOrder(order);
	intrule.Print();
	
	TPZVec<REAL> point2(3,0.L);
	REAL weight2 = 0.L;
	NumInteg intrule2(p);
	intrule2.SetParametersJacobi(0.0L, 0.0L);
	intrule2.SetOrder(order,1);
	intrule2.Print();
	
	unsigned int npoints = intrule.NPoints();
	unsigned int npoints2 = intrule2.NPoints();
	BOOST_CHECK(npoints == npoints2);
	if(npoints != npoints2) return;
	unsigned int it;
	
	// Comparating the points and weights of the two rules
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
		intrule2.Point(it,point2,weight2);
		BOOST_CHECK(fabsl(point[0]-point2[0]) < tol);
		BOOST_CHECK(fabsl(weight-weight2) < tol);
	}
}

template<class T>
void TestingCubatureRuleAllOrders(int type,std::ifstream &olddata) {
	T rule(2);
	int maxord = rule.GetMaxOrder();
	for(int order=0;order<maxord+1;order++) {
		TestingNumericIntegrationRule<T>(expo,order,type,olddata);   // OK
	}
}	
/**
 * For first Function. Using Mathematica, I had the following values:
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

	// Testing over GaussLegendre, GaussLobatto and GaussJacobi rules and over all order < 13
	for(int type=0;type<NTypes;type++) {
		TestingCubatureRuleAllOrders<TPZInt1d>(type,olddata);
		olddata.close();
		olddata.open(filename.c_str());
	}
}

BOOST_AUTO_TEST_CASE(numinteg2DQ_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult2DQ.txt";
	std::ifstream olddata(filename.c_str());

	// Quadrilateral parametric space
	// Testing over GaussLegendre, GaussLobatto and GaussJacobi rules and over all order < 13
	for(int type=0;type<NTypes;type++) {
		TestingCubatureRuleAllOrders<TPZIntQuad>(type,olddata);
		olddata.close();
		olddata.open(filename.c_str());
	}
}

BOOST_AUTO_TEST_CASE(numinteg2DT_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult2DT.txt";
	std::ifstream olddata(filename.c_str());
	
//	int NDigits = NDigitsPrec;
//	NDigitsPrec = 11;
	TestingCubatureRuleAllOrders<TPZIntTriang>(0,olddata);
//	NDigitsPrec = NDigits;
}

BOOST_AUTO_TEST_CASE(numinteg3DC_tests) {

	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3DC.txt";
	std::ifstream olddata(filename.c_str());

	// Cube
	for(int type=0;type<NTypes;type++) {
		TestingCubatureRuleAllOrders<TPZIntCube3D>(type,olddata);
		olddata.close();
		olddata.open(filename.c_str());
	}
}

BOOST_AUTO_TEST_CASE(numinteg3DT_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3DT.txt";
	std::ifstream olddata(filename.c_str());

//	int NDigits = NDigitsPrec;
//	NDigitsPrec = 10;
	TestingCubatureRuleAllOrders<TPZIntTetra3D>(0,olddata);
//	NDigitsPrec = NDigits;
}

BOOST_AUTO_TEST_CASE(numinteg3DPy_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3DPy.txt";
	std::ifstream olddata(filename.c_str());

//	int NDigits = NDigitsPrec;
//	NDigitsPrec = 7;
	TestingCubatureRuleAllOrders<TPZIntPyram3D>(0,olddata);
//	NDigitsPrec = NDigits;
}

BOOST_AUTO_TEST_CASE(numinteg3DPr_tests) {
	
	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
	filename += "FirstResult3DPr.txt";
	std::ifstream olddata(filename.c_str());

//	int NDigits = NDigitsPrec;
	TestingCubatureRuleAllOrders<TPZIntPrism3D>(0,olddata);
//	NDigitsPrec = NDigits;
}

BOOST_AUTO_TEST_SUITE_END()

/**
 * @brief Suite for test of numeric integration of the PolinomialN with lower order and higher order.
 * For lower order the test check must to be right. For higher order the test must to failed.
 */

int NN = 6;
REAL a = -0.5L;
REAL b = 2.1L;

template <class NumInteg>
void TestingNumericLowerIntegrationRule(int order,int type,REAL LimInf,REAL LimSup) {
	// Variables to computing numerical integration
	TPZVec<REAL> point(3,0.L);
	REAL pointx;
	REAL jacob = 0.5L*(b-a);
	REAL weight = 0.0L;
	REAL integral = 0.0L;
	TPZVec<int> orderVec(3,order);
	
	// Integration rule
	NumInteg intrule(order);
	intrule.SetOrder(orderVec,type);
	
	unsigned int npoints = intrule.NPoints();
	unsigned int it;
	
	// Computing the definite integral for Funcao on parametric element
	cout << "Testing PolinomialN : N = " << NN << " \t order = " << order << " \ttype = " << type << " \t NPoints = " << npoints << "  ";
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
		pointx = LinearTransformFromMaster(point[0],LimInf,LimSup);   // x=T(w)   where w pertence [-1.,1.]
		integral +=  weight * PolinomialN(NN,pointx) * jacob;
	}
	
	// Making tol compatible with the wished significant digits
	long double tol = 1.0L;
	REAL integ = IntegralOfPolinomialN(NN,a,b);
	REAL result;
	for(it=0;it<NDigitsPrec;it++)
		tol *= 0.1L;
	result = fabsl(integral-integ);
	if(order <= NN)
		cout << "Difference: " << result << ".\n";
	else {
		cout << std::endl;
		BOOST_CHECK_MESSAGE(result < tol,"Failed Test of Integral with higher order than polinomial: " << result << ".\n");
	}
}

BOOST_AUTO_TEST_SUITE(numint_test2)

BOOST_AUTO_TEST_CASE(LowerOrderTest) {
	NN = 3;
	int NDigits = NDigitsPrec;
	TPZInt1d rule(1);
	int maxord = rule.GetMaxOrder();
	while(NN<(2*maxord)) {
		NN++;
		for(int type=0;type<NTypes;type++) {
			for(int order=0;order<NN;order++) {
				if(order > 10)
					NDigitsPrec = 6;
				TestingNumericLowerIntegrationRule<TPZInt1d>(order+1,type,a,b);
			}
		}
	}
	NDigitsPrec = NDigits;
}

BOOST_AUTO_TEST_CASE(HigherOrderTest) {
	NN = 3;
	int NDigits = NDigitsPrec;
	TPZInt1d rule;
	int maxord = rule.GetMaxOrder();
	while(NN<(2*maxord) && NN<25) {
		NN++;
		for(int type=0;type<NTypes;type++) {
			for(int order=(NN+1);order<(NN+maxord);order++) {
				if(order > 10)
					NDigitsPrec = 6;
				TestingNumericLowerIntegrationRule<TPZInt1d>(order,type,a,b);
			}
		}
	}
	NDigitsPrec = NDigits;
}

BOOST_AUTO_TEST_SUITE_END()

#endif

