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

#include "pzquad.h"

// Using Unit Test of the Boost Library
#ifdef PZ_USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz numericintegration geomtests

#include "boost/test/unit_test.hpp"
#include "boost/test/tools/output_test_stream.hpp"

#endif

std::string dirname = PZSOURCEDIR;
unsigned int NDigitsPrec = 13;
int NTypes = 2;
// For cubature rules
// Conclusion:	Use 4 to run with REAL = float
//				Use 13 to run with REAL = double
//				Use 15 to run with REAL = long double

/** 
 * @name Testing a numeric integration rule of order p for function with argument and coefficients
 * depending on p.
 * @{ 
 */

REAL power(int s, REAL x)
{
    REAL fun = 1.0; for (int is = 0; is < s; is++) { fun *= x;}
    return fun;
}

/**
 * @brief Define a function to integrate, please render it using latex.
 \f$ f\left(\xi,\eta,\zeta\right)=\overset{p}{\underset{i=0}{\sum}}\overset{p}{\underset{j=0}{\sum}}\overset{p}{\underset{j=0}{\sum}}\left(i+1\right)\left(\xi+\frac{1}{10}\right)^{i}\left(j+1\right)\left(\eta+\frac{1}{10}\right)^{j}\left(k+1\right)\left(\zeta+\frac{1}{10}\right)^{k}\;\; with\;\; i+j+k\leq p 
 \f$
*/

REAL Function(TPZManVector<REAL,3> &point, int p, int dim) {
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

/* @} */

/** 
 * @name Testing a numeric integration rule of order p for polinomial of order N
 * @{ 
 */

/** 
 * @brief Compute the value of polinomial \f$ P(x) = (N+1)*x^N + N^2 \f$. \n
 * The integral over \f$ [a,b] \f$ is \f$ \int[P(x)dx]\sub{a} = (b^{N+1} - a^{N+1}) + N^2 (b-a) \f$
 */
REAL PolinomialN(int N,REAL x)
{
	REAL valuex = 1.0L;
    for(int i=0;i<N;i++){ valuex *= x; }
	return (((REAL)(N+1))*valuex + ((REAL)(N*N)));
}

/**
 * @brief Compute the integral value of PolinomialN over the interval \f$[LimInf=a, LimSup=b]\f$ 
 * @param N Degree of the polinomial
 */
REAL IntegralOfPolinomialN(int N,REAL LimInf, REAL LimSup)
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


#ifdef PZ_USING_BOOST

/**
 * @brief Tests the numerical integration rules implemented in Neopz
 * @param order Order of the numeric integration
 * @note NumInteg Respresents the Numeric Integration Rule.
 */


template <class NumInteg>
void TestingNumericIntegrationRule(int p,int type, boost::test_tools::output_test_stream &out) {
	// Variables to computing numerical integration
	TPZManVector<REAL,3> point(3,0.L);
	REAL weight = 0.L;
	REAL NeopzIntegral = 0.L;
	TPZManVector<int,20> order(3,p);

	// Variables to put numerical value as wish format
	char text[64];
	char format[32];
	memset(text,0,strlen(text));
	memset(format,0,strlen(format));
	
	// Creating the integration rule
	NumInteg IntegrationRule(p);
	IntegrationRule.SetType(type);
	IntegrationRule.SetOrder(order);
    int dimension = IntegrationRule.Dimension();

	int npoints = IntegrationRule.NPoints();
    
	// Integrates Fucntion on parametric element space
	for (int it=0;it<npoints;it++) {
		IntegrationRule.Point(it,point,weight);
		NeopzIntegral += weight * Function(point,p,dimension);
	}

	// Formating the integral value obtained
    if(NeopzIntegral < 10.) { sprintf(format,"%%.%dLf",NDigitsPrec);}
    else if(NeopzIntegral < 100.) { sprintf(format,"%%.%dLf",NDigitsPrec-1); }
    else if(NeopzIntegral < 1000.) { sprintf(format,"%%.%dLf",NDigitsPrec-2); }
    else if(NeopzIntegral < 10000.) { sprintf(format,"%%.%dLf",NDigitsPrec-3); }
    else { sprintf(format,"%%.%dLf",NDigitsPrec-4); }
	sprintf(text,format,NeopzIntegral);
	
	// Stores the obtained integral value into the boost::output to compare with match_pattern() i.e. the integral values from mathematica notebook NumInteBeingChecked.nb
    out << text << "\n";
	BOOST_CHECK_MESSAGE( out.match_pattern() , "\nIntegration: Dimension = " << dimension << "\t Order = " << p << "\t Number of Points = " << npoints << "\t Value = " << text << "\n");
    
}
template <class NumInteg>
void TestingNumericIntegrationRule(int p,int type,std::ifstream &input) {
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
	long double tol = 1.L;
	input >> MathematicaIntegral;
#ifdef REALfloat
    NDigitsPrec = 6;
#endif
	// Making tol compatible with the wished significant digits
    for(unsigned int it=0; it < NDigitsPrec; it++){ tol *= 0.1L; }
	
    if(MathematicaIntegral > 10.0) { tol *= 10.L; }
    if(MathematicaIntegral > 100.0) { tol *= 10.L; }
    if(MathematicaIntegral > 1000.0) { tol *= 10.L; }
    REAL result;
#ifndef REALfloat
    result = fabs(NeopzIntegral-MathematicaIntegral);
#else
    result = std::abs(NeopzIntegral-MathematicaIntegral);
#endif
	// If the boolean expresion returns false, then the message will be displayed.
	BOOST_CHECK_MESSAGE(result < tol , "\nIntegration: Dim = " << dimension << "\t Order = " << p << "\t NPoints = " << npoints << "\t Value = " << NeopzIntegral << " difference = " << result << "\n" );
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
	for(int order=0;order <= maxord ;order++) {
		TestingNumericIntegrationRule<T>(order,type,olddata);   // OK
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


BOOST_AUTO_TEST_SUITE(geomtests)


BOOST_AUTO_TEST_CASE(geom_integration_tests) {

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





BOOST_AUTO_TEST_SUITE_END()


#endif

