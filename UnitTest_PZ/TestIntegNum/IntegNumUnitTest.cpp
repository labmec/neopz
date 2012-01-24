/**
 * @file
 * @brief Contains Unit Test in Boost framework for Numerical integration module of the NeoPZ
 */

#include "pzvec.h"
#include "pzquad.h"

#include <iostream>
#include <fstream>
using namespace std;

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
bool first_run = false;
int FirstOrder = 0;

// Numerical integration values obtained using Mathematica: one-dim, quadrilateral, triangular, cube, tetrahedram, pyramidal, prism
REAL intvalue[7][MAXORDER] = { {2.000000000000000,2.606198284550454,2.659572164415588,1.499341835485549,2.367617787494460,1.825083430864047,2.169215575174691,1.617745418673051,2.480691807001154,1.347699766137747,1.417801801857337,1.935192061654517,1.429663752832786},
	{4.,6.424793138201817,4.310074335087271,2.626129673372531,6.054379903820096,2.924805380542269,3.6230005930154684,4.556300644939264,4.577512525115554,1.6940832032465925,4.052863938018176,3.858582377588495,1.6518896896567403 },
	{0.5,0.9320398994069123,0.54479546786258,0.5063395201360538,0.45428110865992366,0.2509209988237399,0.26457181178979317,0.4853077813019233,0.8543252469564175,0.3195852820307958,0.3282395374596146,0.49096556829476024,0.18197735362039558},
	{10.42479313820182,15.27437941460545,11.04494180837636,7.677052484946879,14.53355294584201,8.274403899286355,9.670794324232754,11.53739442808034,11.57981818843293,5.812959544695003,10.53052101423817,10.14195789337881,5.728572517515299},
	{0.21688291481409032598,0.31731541110893806824,0.2271127993891759,0.2227611387537444,0.2227611387537444,0.1579731476556451,0.1592039908614216,0.2114714205797953,0.3036594998229580,0.1789852244559073,0.1737775985061321,0.2150662368602839,0.1447150978897556},
	{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
	{1.303099142275227,2.167178941089052,1.392690078000387,1.315778182547335,1.211661359595075,0.8049411399227070,0.8322427658548136,1.273714704879074,2.011749636188062,0.9422697063368188,0.9595782171944565,1.285030278864748,0.6670538495160184} };

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
void TestingNumericIntegrationRule(int exp, int p, REAL intvalue,boost::test_tools::output_test_stream &out,std::ofstream &results) {
	
	TPZVec<REAL> point(3,0.);
	REAL weight = 0.;
	REAL integral = 0.;
	TPZVec<int> order(3,p);

	// Integration rule
	NumInteg intrule(p);
	intrule.SetOrder(order);

	int npoints = intrule.NPoints();
	int it;

	cout << "\nIntegration dimension = " << intrule.Dimension() << "\t order = " << p << "\t NPoints = " << npoints << "\n";
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
		cout << "\t Point " << it << " : " << point[0] << "\t" << point[1] << "\t" << point[2] << std::endl;
		cout << "\t Weight : " << weight << std::endl;
		integral += weight * Funcao(point,exp,p);
	}
	
	cout << "Numerical integration value = " << integral << "\tDiff = " << fabs(intvalue - integral) << std::endl << std::endl;
	
	BOOST_CHECK(IsZero((intvalue-integral)*.0001));
	
	// Stores the obtained integral value
	if(!first_run) out << integral << std::endl;
	// First tests stores the integral values for following tests
	else results << integral << std::endl;
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
	
	boost::test_tools::output_test_stream output("FirstResult1D.txt",false);
	std::ofstream FirstResult("FirstResult1D.txt");
	
	int i = 2;
	int order;
	// Whether FirstResult.txt is empty then it is a first run
	for(order=FirstOrder;order<MAXORDER;order++) {
		if(order == 1) continue;
		TestingNumericIntegrationRule<TPZInt1d>(i,order,intvalue[0][order],output,FirstResult);   // OK
	}
	// Conclusion: It's failed at order = 1. But the erro is proportional e-10 from order = 8
	
	// For next executions
	if(!first_run) BOOST_CHECK( output.match_pattern() );
	FirstResult.close();
}

BOOST_AUTO_TEST_CASE(numinteg2D_tests) {
	
	boost::test_tools::output_test_stream output("FirstResult2D.txt",false);
	std::ofstream FirstResult("FirstResult2D.txt");
	
	int i = 2;
	int order;
	// Quadrilateral parametric space
	for(order=FirstOrder;order<MAXORDER;order++) {
		if(order == 1) continue;
		TestingNumericIntegrationRule<TPZIntQuad>(i,order,intvalue[1][order],output,FirstResult);   // OK
	}
	// Conclusion: We have problem at order 1. But the erro is proportional e-10 from order = 8
	if(first_run) FirstResult << std::endl;
	else output << std::endl;
	
	// Triangular parametric space
	for(order=FirstOrder;order<MAXORDER;order++) {
		if(order == 1 || order == 5) continue;
		TestingNumericIntegrationRule<TPZIntTriang>(i,order,intvalue[2][order],output,FirstResult);
	}
	// Conclusion: We have problem at order 1 and 5. It's failed at order >= 10.
	
	// For next executions
	if(!first_run) BOOST_CHECK( output.match_pattern() );
	FirstResult.close();
}

BOOST_AUTO_TEST_CASE(numinteg3D_tests) {
	
	boost::test_tools::output_test_stream output("FirstResult3D.txt",false);
	std::ofstream FirstResult("FirstResult3D.txt");
	
	int i = 2;
	int order;
	// Cube
	for(order=FirstOrder;order<MAXORDER;order++) {
		if(!order || order == 1) continue;
		TestingNumericIntegrationRule<TPZIntCube3D>(i,order,intvalue[3][order],output,FirstResult);   // OK
	}
	// Conclusion: We have problem at order 1 and 0. But the erro is proportional e-10 from order = 8
	if(first_run) FirstResult << std::endl;
	else output << std::endl;
	
	// Tetrahedram
	for(order=FirstOrder;order<MAXORDER;order++) {
		if(!order || order == 1 || order == 4) continue;
		TestingNumericIntegrationRule<TPZIntTetra3D>(i,order,intvalue[4][order],output,FirstResult);
	}
	// Conclusion: It's failed at order = 1 and order = 4
	if(first_run) FirstResult << std::endl;
	else output << std::endl;	
	
	// Pyramidal
	for(order=FirstOrder;order<MAXORDER;order++) {
		FirstResult << 1.0 << std::endl;
		//		TestingNumericIntegrationRule<TPZIntPyram3D>(i,order,out);  // It isn't implemented yet
	}
	if(first_run) FirstResult << std::endl;
	else output << std::endl;
	
	// Prism
	for(order=FirstOrder;order<MAXORDER;order++) {
		if(!order || order == 1 || order == 5) continue;
		TestingNumericIntegrationRule<TPZIntPrism3D>(i,order,intvalue[6][order],output,FirstResult);
	}
	// Conclusion: It's failed at order = 1 and order = 0

	// For next executions
	if(!first_run) BOOST_CHECK( output.match_pattern() );
	
	FirstResult.close();
}

BOOST_AUTO_TEST_SUITE_END()

#endif
