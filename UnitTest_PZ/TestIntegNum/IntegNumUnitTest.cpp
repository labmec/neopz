/**
 * @file IntegNumUnitTest.cpp
 * @brief Define a Unit Test using Boost for Numerical integration of the NeoPZ
 *
 */


#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include "pzvec.h"
#include "pzmanvector.h"
#include "pzquad.h"
#include "TPZCurve.h"
#include <catch2/catch.hpp>

static unsigned int NDigitsPrec = 13;
static int NTypes = 2;
// For cubature rules
// Conclusion:	Use 4 to run with REAL = float
//				Use 13 to run with REAL = double
//				Use 15 to run with REAL = long double

/** 
 * @name Testing a numeric integration rule of order p for function with argument and coefficients
 * depending on p.
 * @{ 
 */

static REAL power(int s, REAL x)
{
    REAL fun = 1.0; for (int is = 0; is < s; is++) { fun *= x;}
    return fun;
}

/**
 * @brief Define a function to integrate, please render it using latex.
 \f$ f\left(\xi,\eta,\zeta\right)=\overset{p}{\underset{i=0}{\sum}}\overset{p}{\underset{j=0}{\sum}}\overset{p}{\underset{j=0}{\sum}}\left(i+1\right)\left(\xi+\frac{1}{10}\right)^{i}\left(j+1\right)\left(\eta+\frac{1}{10}\right)^{j}\left(k+1\right)\left(\zeta+\frac{1}{10}\right)^{k}\;\; with\;\; i+j+k\leq p 
 \f$
*/

static REAL Function(TPZManVector<REAL,3> &point, int p, int dim);
/* @} */

/** 
 * @name Testing a numeric integration rule of order p for polinomial of order N
 * @{ 
 */

/** 
 * @brief Compute the value of polinomial \f$ P(x) = (N+1)*x^N + N^2 \f$. \n
 * The integral over \f$ [a,b] \f$ is \f$ \int[P(x)dx]\sub{a} = (b^{N+1} - a^{N+1}) + N^2 (b-a) \f$
 */
static REAL PolinomialN(int N,REAL x)
{
	REAL valuex = 1.0L;
    for(int i=0;i<N;i++){ valuex *= x; }
	return (((REAL)(N+1))*valuex + ((REAL)(N*N)));
}

/**
 * @brief Compute the integral value of PolinomialN over the interval \f$[LimInf=a, LimSup=b]\f$ 
 * @param N Degree of the polinomial
 */
static REAL IntegralOfPolinomialN(int N,REAL LimInf, REAL LimSup)
{
	REAL A = 1.0L, B = 1.0L;
    for(int i=0;i<(N+1);i++) { A *= LimInf; B *= LimSup; }
	return ((B-A)+(((REAL)(N*N))*(LimSup-LimInf)));
}

/**
 * @brief Linear transformation to transform a interval [-1.,1.] to [a,b] 
 * @param point Point within the master element [-1.,1.]
 * @param LimInf Initial point of the interval image of the master
 * @param LimSup Final point of the interval image of the master
 */
static REAL LinearTransformFromMaster(REAL point,REAL LimInf, REAL LimSup) {
	return (0.5L * ((LimSup - LimInf)*point + (LimSup + LimInf)));
}

/**
 * @brief Linear transformation to transform a interval [-1.,1.] to [a,b] 
 * @param point Point within the master element [-1.,1.]
 * @param LimInf Initial point of the interval image of the master
 * @param LimSup Final point of the interval image of the master
 */
static REAL DifferenceOfLinearTransformFromMaster(REAL point,REAL LimInf, REAL LimSup) {
	return (0.5L * (LimSup - LimInf));
}

/** @} */


/**
 * @brief Tests the numerical integration rules implemented in Neopz
 * @param order Order of the numeric integration
 * @note NumInteg Respresents the Numeric Integration Rule.
 */
template <class NumInteg>
void TestingNumericIntegrationRule(int p,int type,std::ifstream &input);

template<class NumInteg>
void ComparePointsWithNewRule(int p);

template<class T>
void TestingCubatureRuleAllOrders(int type,std::ifstream &olddata) {
	T rule(2);
    int maxord = rule.GetMaxOrder();
	for(int order=0;order <= maxord ;order++) {
		TestingNumericIntegrationRule<T>(order,type,olddata);   // OK
	}
}

void ComputeError(REAL alpha, TPZManVector<REAL,3> &coordinate,TPZGeoEl * GeometricEl, REAL &error);

void TaylorCheck(TPZManVector<REAL,3> &coordinate,TPZGeoEl * GeometricEl);

void IntegrateCurve(TPZCurve &curve);

TEST_CASE("numinteg1D_tests","[numinteg_tests]") {

	// File with integration values calculated previously
	std::string filename = "Line.txt";
	std::ifstream MathematicaData(filename.c_str());

	// Testing over Linbo Zhang rule and over all order < 13
	for(int type=0;type<NTypes;type++) {
		TestingCubatureRuleAllOrders<TPZInt1d>(type,MathematicaData);
		MathematicaData.close();
		MathematicaData.open(filename.c_str());
	}
}

//TEST_CASE("numinteg2DQ_tests","[numinteg_tests]") {
//	
//	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
//	filename += "Quadrilateral.txt";
//	std::ifstream MathematicaData(filename.c_str());
//
//	// Testing over GaussLegendre, GaussLobatto and GaussJacobi rules and over all order < 13
//	for(int type=0;type<NTypes;type++) {
//		TestingCubatureRuleAllOrders<TPZIntQuad>(type,MathematicaData);
//		MathematicaData.close();
//		MathematicaData.open(filename.c_str());
//	}
//}

TEST_CASE("numinteg2DT_tests","[numinteg_tests]") {
	
	std::string filename = "Triangle.txt";
	std::ifstream MathematicaData(filename.c_str());
    // Testing over GaussLegendre, GaussLobatto and GaussJacobi rules and over all order < 13
    const int type = 0;// triangle has only one integration rule
    TestingCubatureRuleAllOrders<TPZIntTriang>(type,MathematicaData);
    MathematicaData.close();
    MathematicaData.open(filename.c_str());

}

//TEST_CASE("numinteg3DC_tests","[numinteg_tests]") {
//
//	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
//	filename += "Cube.txt";
//	std::ifstream MathematicaData(filename.c_str());
//    
//	for(int type=0;type<NTypes;type++) {
//		TestingCubatureRuleAllOrders<TPZIntCube3D>(type,MathematicaData);
//		MathematicaData.close();
//		MathematicaData.open(filename.c_str());
//	}
//}
// 
//TEST_CASE("numinteg3DT_tests","[numinteg_tests]") {
//	
//	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
//	filename += "Tetrahedral.txt";
//	std::ifstream MathematicaData(filename.c_str());
//	TestingCubatureRuleAllOrders<TPZIntTetra3D>(0,MathematicaData);
//    MathematicaData.close();
//}
//
//TEST_CASE("numinteg3DPy_tests","[numinteg_tests]") {
//	
//	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
//	filename += "Pyramid.txt";
//	std::ifstream MathematicaData(filename.c_str());
//	TestingCubatureRuleAllOrders<TPZIntPyram3D>(0,MathematicaData);
//    MathematicaData.close();
//}
//
//TEST_CASE("numinteg3DPr_tests","[numinteg_tests]") {
//	
//	std::string filename = dirname + "/UnitTest_PZ/TestIntegNum/";
//	filename += "Prism.txt";
//	std::ifstream MathematicaData(filename.c_str());
//	TestingCubatureRuleAllOrders<TPZIntPrism3D>(0,MathematicaData);
//    MathematicaData.close();
//}



///**
// * @brief Suite for test of numeric integration of the PolinomialN with lower order and higher order.
// * For lower order the test check must to be right. For higher order the test must to failed.
// */
//
//
//int NN = 10;
//REAL a = -0.5L;
//REAL b = 2.1L;
//
//template <class NumInteg>
//void TestingNumericLowerIntegrationRule(int order,int type,REAL LimInf,REAL LimSup) {
//	// Variables to computing numerical integration
//	TPZManVector<REAL,3> point(3,0.L);
//	REAL pointx;
//	REAL jacob = 0.5L*(b-a);
//	REAL weight = 0.0L;
//	REAL NeopzIntegral = 0.0L;
//	TPZManVector<int,20> orderVec(3,order);
//	
//	// Creating the integration rule
//	NumInteg IntegrationRule(order);
//	IntegrationRule.SetOrder(orderVec,type);
//	
//	unsigned int npoints = IntegrationRule.NPoints();
//	
//	// Computing the definite integral for Funcao on parametric element
//	for (unsigned int it =0;it<npoints;it++) {
//		IntegrationRule.Point(it,point,weight);
//		pointx = LinearTransformFromMaster(point[0],LimInf,LimSup);
//		NeopzIntegral +=  weight * PolinomialN(order,pointx) * jacob;
//	}
//	
//	// Making tol compatible with the wished significant digits
//	long double tol = 1.0L;
//	REAL Integral = IntegralOfPolinomialN(order,a,b);
//	REAL result;
//    for(unsigned int it=0; it < NDigitsPrec; it++){	tol *= 0.1L; }
//    result = fabsl(NeopzIntegral-Integral);
//	
////    if(order <= NN)
////        std::cout << "Difference: " << result << ".\n";
////	else {
////        std::cout << std::endl;
////        std::cout<<"Failed Test of Integral with higher order than polinomial: " << result << ".\n";
////		REQUIRE(result < tol);
////	}
//    if(order >= NN)
//    {
//        std::cout << "Difference: " << result << ".\n";
//        std::cout << std::endl;
//        // If the boolean expresion returns false, then the message will be displayed.
//        std::cout<<"Failed Test of Integral with higher order than polinomial: " << result << ".\n";
//        REQUIRE(result < tol);
//    }
// 
//}
//
//
////TEST_CASE("LowerOrderTest","[numint_test2]") {
////	int NDigits = NDigitsPrec;
////	TPZInt1d rule(1);
////	int maxord = rule.GetMaxOrder();
////	while(NN<(2*maxord)) {
////		NN++;
////		for(int type=0;type<1;type++) {
////			for(int order=0;order<NN;order++) {
////				if(order > 10)
////					NDigitsPrec = 6;
////				TestingNumericLowerIntegrationRule<TPZInt1d>(order+1,type,a,b);
////			}
////		}
////	}
////	NDigitsPrec = NDigits;
////}
////
////TEST_CASE("HigherOrderTest","[numint_test2]") {
////	int NDigits = NDigitsPrec;
////	TPZInt1d rule;
////	int maxord = rule.GetMaxOrder();
////	while(NN<(2*maxord) && NN<25) {
////		NN++;
////		for(int type=0;type<1;type++) {
////			for(int order=(NN+1);order<(NN+maxord);order++) {
////				if(order > 10)
////					NDigitsPrec = 6;
////				TestingNumericLowerIntegrationRule<TPZInt1d>(order,type,a,b);
////			}
////		}
////	}
////	NDigitsPrec = NDigits;
////}
//



TEST_CASE("geom_integration_tests","[geomtests]") {

    int href = 5;
//  Check taylor convergence for all curve elements
    TPZCurve * Curve = new TPZCurve;
    Curve->SetRadius(1.0);

    Curve->MakeRhombus();
    Curve->RefineMe(href);
    Curve->PrintMe();
    IntegrateCurve(*Curve);

    Curve->MakeCircleWave();
    Curve->RefineMe(href);
    Curve->PrintMe();
    IntegrateCurve(*Curve);
    
    Curve->MakeCircleFromArc();
    Curve->RefineMe(href);
    IntegrateCurve(*Curve);

}


static REAL Function(TPZManVector<REAL,3> &point, int p, int dim) {
	REAL functionvalue = (REAL)(0.0L);
	REAL xi = point[0];
	REAL eta = point[1];
	REAL zeta = point[2];
    REAL epsilon = 0.1;
    int p1 = 0;
    int p2 = 0;
    int p3 = 0;
    
    switch (dim) {
        case 1:
        {   p1 = p; }
            break;
        case 2:
        {   p1 = p; p2 = p; }
            break;
        case 3:
        {   p1 = p; p2 = p; p3 = p; }
            break;
        default:
        {   std::cout << "Inconsistency: wrong dimension " << dim << std::endl;
            DebugStop();
        }
            break;
    }
    
    for(int k = 0; k <= p3; k++)
    {
        for(int j = 0; j <= p2; j++)
        {
            for(int i = 0; i <= p1; i++)
            {
                if (i+j+k <= p) {
                    functionvalue += (i+1)*power(i,xi+epsilon) * (j+1)*power(j,eta+epsilon) * (k+1)*power(k,zeta+epsilon);
                }
            }
        }
    }
	return functionvalue;
}


template <class NumInteg>
void TestingNumericIntegrationRule(int p,int type,std::ifstream &input) {
  // Check if file exists
  if (!input.is_open()) {
    std::cout << "Error: Input file not found!\n";
    DebugStop();
  }

	// Variables to computing numerical integration
	TPZManVector<REAL,3> point(3);
	REAL weight = 0.L;
	REAL NeopzIntegral = 0.L;
	TPZManVector<int,20> order(3,p);
		
	// Creating the integration rule
	NumInteg IntegrationRule(p);
	IntegrationRule.SetOrder(order,type);
	
	unsigned int npoints = IntegrationRule.NPoints();
    int dimension = IntegrationRule.Dimension();
	std::string namerule;
	IntegrationRule.Name(namerule);
    
	
	// Integrates Fucntion on parametric element space
    std::cout << " Cubature rule: " << namerule << "  Type " << type << " \tOrder " << p << " \tNumber of Points " << npoints << std::endl;
	for (unsigned int it = 0;it<npoints;it++) {
		IntegrationRule.Point(it,point,weight);
		NeopzIntegral += weight * Function(point,p,dimension);
	}

	// Variables to import data from files generate by notebook NumInteBeingChecked.nb
	long double MathematicaIntegral;
	long double tol = 6.L;
	input >> MathematicaIntegral;
	// Making tol compatible with the wished significant digits
#ifdef REALfloat
    NDigitsPrec = 6;
#endif
    for(unsigned int it=0; it < NDigitsPrec; it++){ tol *= 0.1L; }
	
    if(MathematicaIntegral > 10.0) { tol *= 10.L; }
    if(MathematicaIntegral > 100.0) { tol *= 10.L; }
    if(MathematicaIntegral > 1000.0) { tol *= 10.L; }
    REAL result = std::abs(NeopzIntegral-MathematicaIntegral);
	// If the boolean expresion returns false, then the message will be displayed.
    if(!(result<tol)){
      std::cerr << "\nIntegration: Dim = " << dimension << "\t Order = " << p
                << "\t NPoints = " << npoints << "\t Value = " << NeopzIntegral
                << " difference = " << result << "\n";
    }
	REQUIRE(result < tol);
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
	REQUIRE(npoints == npoints2);
	if(npoints != npoints2) return;
	unsigned int it;
	
	// Comparating the points and weights of the two rules
	for (it=0;it<npoints;it++) {
		intrule.Point(it,point,weight);
		intrule2.Point(it,point2,weight2);
		REQUIRE(fabsl(point[0]-point2[0]) < tol);
		REQUIRE(fabsl(weight-weight2) < tol);
	}
}

void ComputeError(REAL alpha, TPZManVector<REAL,3> &coordinate,TPZGeoEl * GeometricEl, REAL &error){

  int dimension = GeometricEl->Dimension();
  TPZFMatrix<REAL> GradofX;
  TPZFMatrix<REAL> GradofXAlpha;
  TPZFMatrix<REAL> Error;
  TPZFMatrix<REAL> X(3,1,0.0);
  TPZFMatrix<REAL> XAlpha(3,1,0.0);
  TPZFMatrix<REAL> Alpha(dimension,1,alpha);

  TPZFMatrix<REAL> jac;
  TPZFMatrix<REAL> axes;
  TPZFMatrix<REAL> axesT;
  REAL detjac;
  TPZFMatrix<REAL> jacinv;

  GeometricEl->Jacobian(coordinate,jac,axes,detjac,jacinv);
  axes.Transpose(&axesT);
  axesT.Multiply(jac,GradofX);

  TPZManVector<REAL,3> coordinateAlpha(coordinate);
  TPZManVector<REAL,3> result(3,0.0);
  TPZManVector<REAL,3> resultAlpha(3,0.0);

  coordinateAlpha[0] += alpha;
  coordinateAlpha[1] += alpha;
  coordinateAlpha[2] += alpha;

  GeometricEl->X(coordinate,result);
  GeometricEl->X(coordinateAlpha,resultAlpha);

  X(0,0) = result[0];
  X(1,0) = result[1];
  X(2,0) = result[2];

  XAlpha(0,0) = resultAlpha[0];
  XAlpha(1,0) = resultAlpha[1];
  XAlpha(2,0) = resultAlpha[2];

  GradofX.Multiply(Alpha,GradofXAlpha);
  Error = XAlpha - (X + GradofXAlpha);
  error = Norm(Error);

}

void TaylorCheck(TPZManVector<REAL,3> &coordinate,TPZGeoEl * GeometricEl)
{
  
  REAL alpha1 = 0.01;
  REAL alpha2 = 0.1;
  REAL error1 = 0.0;
  REAL error2 = 0.0;
  REAL epsilon = 1.0e-4;
  
  ComputeError(alpha1,coordinate,GeometricEl,error1);
  ComputeError(alpha2,coordinate,GeometricEl,error2);

  REAL m = 0.0;
  if(GeometricEl->IsLinearMapping()){
    m = 2.0;
  }
  else{
      m = (log(error1)-log(error2))/(log(alpha1)-log(alpha2));
  }

  if((2.0 - epsilon) > m && m > (2.0 + epsilon) ){
      std::cout << "Error: Wrong in expected approximation Rate m = " << m << std::endl;
  }

  
}

void IntegrateCurve(TPZCurve &curve)
{
  int IntegrationOrder = 0;
  int type = 0;  
  TPZVec<int> order(3,IntegrationOrder);
  
  TPZManVector<REAL,3> point(3,0.L);
  TPZFMatrix<REAL> jac;
  TPZFMatrix<REAL> axes;
  REAL detjac;
  TPZFMatrix<REAL> jacinv;

  TPZGeoMesh * gmesh = curve.GetGeometry();
  int64_t NumberofElements	 = gmesh->NElements();
  STATE Length = 0.0;
  
  for(int iel = 0; iel < NumberofElements; iel++)
  {
    TPZGeoEl * gel = gmesh->Element(iel);
    if(!gel) { continue; }
    if(gel->HasSubElement()) {continue;}
    int geldim = gel->Dimension();

    if(geldim != 1){
      std::cout << "Only one-dimensional elements are integrated." << std::endl;
      DebugStop();      
    }
    
    // Creating the integration rule
    TPZInt1d IntegrationRule(IntegrationOrder);
    IntegrationRule.SetType(type,IntegrationOrder);
    int npoints = IntegrationRule.NPoints();
    REAL weight = 0.0;

    // Integrates Fucntion on parametric element space
    for (int it=0;it<npoints;it++) {
      IntegrationRule.Point(it,point,weight);
      gel->Jacobian(point,jac,axes,detjac,jacinv);
      TaylorCheck(point,gel);
      Length += weight * detjac ;
    }
    
  }

  std::cout << "Curve length = " << Length << std:: endl;
  
}