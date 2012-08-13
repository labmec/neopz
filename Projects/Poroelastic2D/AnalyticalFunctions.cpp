/*
 *  AnalyticalFunctions.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "AnalyticalFunctions.h"
#ifdef USING_BOOST
	#include <boost/math/special_functions/erf.hpp> //Required for erfc function on windows
#endif

//	Internally Defined Constants
static const long double epsilon = 10.0 * LDBL_EPSILON;

using namespace std;


	// Right handside term of our Linear PDE
	void ForcingTimeDependFunction(TPZVec<REAL> &ptx, REAL TimeValue,int WhichStateVariable,double &StateVariable) 
	{
		
		// Define the relations for each variable in the right hand side of the StateVariable at the current PDE.
		
		REAL x = ptx[0];
		REAL y = ptx[1];
		REAL z = 0.0;
		
		REAL hour = 3600;
		REAL day = 86400;
		REAL month = 30*day;
		REAL year = 365*day;
		REAL delta = 99.9999*hour;
		REAL MaxTime = 100.0*hour;
		
		
		switch (WhichStateVariable) 
		{
			case 0:
			{
				//	Ux
				StateVariable = 0.0;
				break;
			}
			case 1:
			{
				//	Uy
				StateVariable = 0.0;
				break;			
			}
			case 2:
			{
				//	Pressure
				//			if ((abs(x-347.15922101486848) < 1.0e-4) && (abs(y-347.15922101486848) < 1.0e-4)) 
				//			{
				//			StateVariable = -0.25*5.0e-5;
				//			}
				//			else 
				//			{
				StateVariable = 0.0;
				//			}
				break;
			}
			default:
				break;
		}
	}


	void ExactSolution1D(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
		//REAL x = ptx[0];
		REAL y = ptx[1];
		
		REAL pini = 1000.;
		REAL lamb = 8333.33;
		REAL mi = 12500.0;
		REAL visc =0.001; 
		REAL perm =  1.e-10;
		REAL H=1.;
		REAL segtime = 0.0;						//	[s]
		segtime = timestep;		
		
		
		int in;
		REAL pD = 0.0, uD = 0.0, sigD = 0.0;
		REAL PI = atan(1.)*4.;
		
		sol[0]=0.0;
		sol[1]=0.0;
		sol[2]=0.0;
		
		REAL tD = (lamb+2.*mi)*perm*segtime/(visc*H);
		REAL xD = abs(y-1.0)/H;
		REAL M = 0.0;
		for (in =0; in<1000; in++) {
			
			M = PI*(2.*in+1.)/2.;
			pD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
			uD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
			sigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
		}
		
		sol[0] = pD*pini;
		sol[1] = (-1.+ sigD)*pini;
		sol[2] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
	}


	void ExactSolution2DLinesource(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
		
		// Defition of varibles
		REAL x = ptx[0];//-500000.0;
		REAL y = ptx[1];//-500000.0;
		REAL z = 0.0;
		REAL r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
		if ( r == 0.0) {
			x=1.0e-20;
			y=1.0e-20;
			r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
		}
		
		
		// Definitions of Parameters
		REAL lamb = 8.3e9;						//	[Pa]
		REAL lambu = 13.4e9;						//	[Pa]
		REAL alpha = 0.7;						//	[-]
		REAL G= 5.5e9;							//	[Pa]
		REAL rhof = 1000.0;						//	[kg/m3]	
		REAL Ssig = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
		REAL segtime = 0.0;						//	[s]
		REAL qMod = -20.0;						//	[kg/s]
		REAL cMod = 0.083;						//	[m2/s]
		REAL kdif = cMod*Ssig;
		REAL PI = atan(1.)*4.;
		segtime = timestep;
		
		if (segtime == 0.0) {
			segtime = 1.0e-20;
		}
		
		
		sol[0]=0.;		// Pressure
		sol[1]=0.;		// Ux
		sol[2]=0.;		// Uy
		flux(0,0)=0.0; // SigX
		flux(1,0)=0.0; // SigY
		flux(2,0)=0.0; // TauXY
		
		
		REAL Pressure = 0.0;					//	[Pa]
		REAL Ux = 0.0;							//	[m]
		REAL Uy = 0.0;							//	[m]
		REAL Sigx = 0.0;						//	[Pa]
		REAL Sigy = 0.0;						//	[Pa]
		REAL Tauxy = 0.0;						//	[Pa]
		
		REAL Zz = (pow(r, 2)/(4.0*cMod*segtime));
		REAL Ei = Exponential_Integral_Ei(-Zz);
		REAL Den = (8*PI*rhof*kdif*(lamb+2.0*G));
		REAL Nem = qMod*alpha*x;
		
		Pressure = (qMod/(4*PI*rhof*kdif))*Exponential_Integral_Ei(-Zz);
		Ux = ((-qMod*alpha*x)/(8*PI*rhof*kdif*(lamb+2.0*G)))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
		Uy = ((-qMod*alpha*y)/(8*PI*rhof*kdif*(lamb+2.0*G)))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
		Sigx = (qMod*alpha*G/(4*PI*rhof*kdif*(2.0*G+lamb)))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(y,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
		Sigy = (qMod*alpha*G/(4*PI*rhof*kdif*(2.0*G+lamb)))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(x,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
		Tauxy = (2.0*qMod*alpha*G*x*y/(4*PI*rhof*kdif*(2.0*G+lamb)*pow(r,2)))*((1/Zz)*(1-exp(-Zz)));
		
		sol[0] = Pressure;
		sol[1] = Ux;
		sol[2] = Uy;
		
		flux(0,0)=Sigx; // SigX
		flux(1,0)=Sigy; // SigY
		flux(2,0)=Tauxy; // TauXY
	}

	void SolucaoExataRosa1D(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
		//REAL x = ptx[0];
		REAL x = ptx[0];//-12500;
		
		
		//	Parameters
		REAL Eyoung = 1.43e10;					//	[Pa]
		REAL poisson = 0.3;						//	[-]
		REAL alpha=0.0;///((1.19667e10)*3);//1.70667e10		1.19667e10					//	[-]
		
		REAL Phi= 1.0;//9.60784e-11;//1.35695e-10				//	[-]
		REAL Ct = 1.0;					//	[kg/m3]
		REAL Se = Ct*Phi;							//	[m2/s]	
		REAL Visc = 1.0;						//	[Pa.s]
		REAL Kmed = 1.0;//7.8784288e-15;//1.109542e-14						//	[m2]
		REAL Qo = 2.0;
		REAL Bo = 1.0;
		REAL PI = atan(1.)*4.;
		REAL segtime = 0.0;					//	[s]	
		
		
		//	REAL Phi= 0.18;//9.60784e-11;//1.35695e-10				//	[-]
		//	REAL Ct = (150.0e-6)*(1/(98066.50));					//	[kg/m3]
		//	REAL Se = Ct*Phi;							//	[m2/s]	
		//	REAL Visc = 0.8*(1.e-3);						//	[Pa.s]
		//	REAL Kmed = 20*(9.86923e-16);//7.8784288e-15;//1.109542e-14						//	[m2]
		//	REAL Qo = 400.0/(86400);
		//	REAL Bo = 1.2;
		//	REAL PI = atan(1.)*4.;
		//	REAL segtime = 0.0;					//	[s]		
		
		segtime = timestep;
		
		if (segtime == 0.0) {
			segtime = 1.0e-12;
		}
		
		x = abs(x);
		sol[0]=0.;
		sol[1]=0.;
		sol[2]=0.;
		
		REAL Eta = (Kmed)/(Phi*Visc*Ct);
		#ifdef USING_BOOST
			REAL Pressure = ((0.5*Qo*Bo*Visc)/(Kmed*1.0))*(sqrt((4*Eta*segtime)/PI)*(exp(-(pow(x,2)/(4*Eta*segtime))))-(x*boost::math::erfc(x/sqrt(4*Eta*segtime))));
		#else
			REAL Pressure = ((0.5*Qo*Bo*Visc)/(Kmed*1.0))*(sqrt((4*Eta*segtime)/PI)*(exp(-(pow(x,2)/(4*Eta*segtime))))-(x*erfc(x/sqrt(4*Eta*segtime))));
		#endif
		
		sol[0] = Pressure;
		//	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
		//	sol[2] = (-1.+ sigD)*pini;
	}

	// Analytical solution for flamant problem Elasticity Theory: Applications and Numerics 
	void ExactSolutionFlamantProblem(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux)
	{
		
		// Defition of variables
		REAL x = ptx[0]-50000.0;
		REAL y = ptx[1]-50000.0;
		REAL z = 0.0;
		REAL Yforce = 1000.0;
		REAL PI = atan(1.)*4.;	
		REAL r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
		if ( r == 0.0) {
			x=1.0e-10;
			y=1.0e-10;
			r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
		}
		
		
		// Definitions of Parameters
		//	REAL lamb = 1.0e9;						//	[Pa]
		//	REAL G= 1.0e9;							//	[Pa]
		//	REAL rhof = 1.0;						//	[kg/m3]	
		//	REAL qMod = -1.0;						//	[kg/s]
		//	REAL cMod = 1.0;						//	[m2/s]
		//	REAL kdif = cMod*Se;
		//	REAL PI = atan(1.)*4.;
		REAL sigXX = -2.0*((Yforce*pow(x,2)*y)/(PI*(pow(r,2))));
		REAL sigYY = -2.0*((Yforce*pow(y,3))/(PI*(pow(r,2))));
		REAL tauXY = -2.0*((Yforce*pow(y,2)*x)/(PI*(pow(r,2))));
		
		sol[0] = sigXX;
		sol[1] = sigYY;
		sol[2] = tauXY;
	}


	void SolucaoExataRosa1DPseudo(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
		//REAL x = ptx[0];
		REAL x = ptx[0];//-12500;
		REAL L = 20000.0;
		//	x = abs(x);
		sol[0]=0.;
		sol[1]=0.;
		sol[2]=0.;
		
		REAL Pressure = 1000 + L*((x/L)-0.5*(pow((x/L),2)));
		
		sol[0] = Pressure;
		//	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
		//	sol[2] = (-1.+ sigD)*pini;
	}


////////////////////////////////////////////////////////////////////////////////
//    Exponential_Integral_Ei                                                 //
//    xExponential_Integral_Ei                                                //
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Ei( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Exponential_Integral_Ei( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Exponential_Integral_Ei( double x )
{
	return (double) xExponential_Integral_Ei( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xExponential_Integral_Ei( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral Ei().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xExponential_Integral_Ei( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xExponential_Integral_Ei( long double x )
{
	if ( x < -5.0L ) return Continued_Fraction_Ei(x);
	if ( x == 0.0L ) return -DBL_MAX;
	if ( x < 6.8L )  return Power_Series_Ei(x);
	if ( x < 50.0L ) return Argument_Addition_Series_Ei(x);
	return Continued_Fraction_Ei(x);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_Ei( long double x )                  //
//                                                                            //
//  Description:                                                              //
//     For x < -5 or x > 50, the continued fraction representation of Ei      //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of Ei(x) is:                          //
//        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_Ei( long double x )
{
	long double Am1 = 1.0L;
	long double A0 = 0.0L;
	long double Bm1 = 0.0L;
	long double B0 = 1.0L;
	long double a = expl(x);
	long double b = -x + 1.0L;
	long double Ap1 = b * A0 + a * Am1;
	long double Bp1 = b * B0 + a * Bm1;
	int j = 1;
	
	a = 1.0L;
	while ( fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1) ) {
		if ( fabsl(Bp1) > 1.0L) {
			Am1 = A0 / Bp1;
			A0 = Ap1 / Bp1;
			Bm1 = B0 / Bp1;
			B0 = 1.0L;
		} else {
			Am1 = A0;
			A0 = Ap1;
			Bm1 = B0;
			B0 = Bp1;
		}
		a = -j * j;
		b += 2.0L;
		Ap1 = b * A0 + a * Am1;
		Bp1 = b * B0 + a * Bm1;
		j += 1;
	}
	return (-Ap1 / Bp1);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_Ei( long double x )                        //
//                                                                            //
//  Description:                                                              //
//     For -5 < x < 6.8, the power series representation for                  //
//     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
//     constant.                                                              //
//     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
//     returned.                                                              //
//                                                                            //
//     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
//        - Sum(1 + 1/2 + ... + 1/j) (-x)^j / j!, where the Sum extends       //
//        from j = 1 to inf.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_Ei( long double x )
{ 
	long double xn = -x;
	long double Sn = -x;
	long double Sm1 = 0.0L;
	long double hsum = 1.0L;
	long double g = 0.5772156649015328606065121L;
	long double y = 1.0L;
	long double factorial = 1.0L;
	
	if ( x == 0.0L ) return (long double) -DBL_MAX;
	
	while ( fabsl(Sn - Sm1) > epsilon * fabsl(Sm1) ) {
		Sm1 = Sn;
		y += 1.0L;
		xn *= (-x);
		factorial *= y;
		hsum += (1.0 / y);
		Sn += hsum * xn / factorial;
	}
	return (g + logl(fabsl(x)) - expl(x) * Sn);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Argument_Addition_Series_Ei(long double x)              //
//                                                                            //
//  Description:                                                              //
//     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
//     Ei.                                                                    //
//                                                                            //
//     The argument addition series for Ei(x) is:                             //
//     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
//     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
//     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
//     from k = 0 to k = j.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////
static long double Argument_Addition_Series_Ei(long double x)
{
	static long double ei[] = {
		1.915047433355013959531e2L,  4.403798995348382689974e2L,
		1.037878290717089587658e3L,  2.492228976241877759138e3L,
		6.071406374098611507965e3L,  1.495953266639752885229e4L,
		3.719768849068903560439e4L,  9.319251363396537129882e4L,
		2.349558524907683035782e5L,  5.955609986708370018502e5L,
		1.516637894042516884433e6L,  3.877904330597443502996e6L,
		9.950907251046844760026e6L,  2.561565266405658882048e7L,
		6.612718635548492136250e7L,  1.711446713003636684975e8L,
		4.439663698302712208698e8L,  1.154115391849182948287e9L,
		3.005950906525548689841e9L,  7.842940991898186370453e9L,
		2.049649711988081236484e10L, 5.364511859231469415605e10L,
		1.405991957584069047340e11L, 3.689732094072741970640e11L,
		9.694555759683939661662e11L, 2.550043566357786926147e12L,
		6.714640184076497558707e12L, 1.769803724411626854310e13L,
		4.669055014466159544500e13L, 1.232852079912097685431e14L,
		3.257988998672263996790e14L, 8.616388199965786544948e14L,
		2.280446200301902595341e15L, 6.039718263611241578359e15L,
		1.600664914324504111070e16L, 4.244796092136850759368e16L,
		1.126348290166966760275e17L, 2.990444718632336675058e17L,
		7.943916035704453771510e17L, 2.111342388647824195000e18L,
		5.614329680810343111535e18L, 1.493630213112993142255e19L,
		3.975442747903744836007e19L, 1.058563689713169096306e20L
	};
	int  k = (int) (x + 0.5);
	int  j = 0;
	long double xx = (long double) k;
	long double dx = x - xx;
	long double xxj = xx;
	long double edx = expl(dx);
	long double Sm = 1.0L;
	long double Sn = (edx - 1.0L) / xxj;
	long double term = DBL_MAX;
	long double factorial = 1.0L;
	long double dxj = 1.0L;
	
	while (fabsl(term) > epsilon * fabsl(Sn) ) {
		j++;
		factorial *= (long double) j;
		xxj *= xx;
		dxj *= (-dx);
		Sm += (dxj / factorial);
		term = ( factorial * (edx * Sm - 1.0L) ) / xxj;
		Sn += term;
	}
	
	return ei[k-7] + Sn * expl(xx); 
}
