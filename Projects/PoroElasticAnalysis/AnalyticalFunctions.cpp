/*
 *  AnalyticalFunctions.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

//#include "specialFunctions.h"

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
		REAL x, y, z;
		x = ptx[0];
		y = ptx[1];
		z = 0.0;
		
		REAL hour = 3600;
		REAL day = 86400;
        REAL month, year, delta, MaxTime;
		month = 30*day;
		year = 365*day;
		delta = 99.9999*hour;
		MaxTime = 100.0*hour;
		
		
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


	void ExactSolutionfiniteColumn1D(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
		//REAL x = ptx[0];
		REAL y = ptx[1];
		
		// Definitions	
		REAL PI = atan(1.)*4.;
        REAL c, cm, SigTop, Gamma, rockrho, gravity, visc, Po;
		c = 1.6;							//	[m2/s]
		cm = 4.9375e-11;					//	[1/Pa]
		SigTop = 1.e8;						//	[Pa]
		Gamma = 0.3333333333;					//	[-]
		rockrho = 2500.0;					//	[kg/m3]
		gravity = 9.81;					//	[m/s2]
		visc =0.001;
		Po = SigTop*Gamma;
		
		REAL segtime = 0.0;						//	[s]
		segtime = timestep;		
		
		
		int in;
		REAL pD = 0.0, uD = 0.0, sigD = 0.0, qDxD = 0.0;
		
		sol[0]=0.0;
		sol[1]=0.0;
		sol[2]=0.0;

		REAL H=1.0;		
		REAL tD = segtime;
		REAL xD = abs(y)/H;
		REAL M = 0.0;
		for (in =0; in<1000; in++) {
			
			M = PI*(2.*in+1.)/2.;
			pD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
			qDxD += -(2.0)*cos(M*xD)*exp(-1.*M*M*tD);
			uD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
			sigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
		}
				
		sol[0] = pD;
		sol[1] = 0.0;
		sol[2] = -(1.- xD - uD);	
				
		flux(0,0)=0.0;
		flux(1,0)=-(1.- sigD);
		flux(2,0)=0.0;
		flux(3,0)=0.0;
		flux(4,0)=-qDxD;
		
	}

void ExactSolutionSemiInfiniteColumn1D(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux)
{
	//REAL x = ptx[0];
	REAL y = ptx[1];
	
	// Definitions
	REAL PI = atan(1.)*4.;
	REAL c, cm, SigTop, Gamma, rockrho, gravity, visc, Po, segtime;
	c = 1.6;							//	[m2/s]
	cm = 4.9375e-11;					//	[1/Pa]
	SigTop = 1.e8;						//	[Pa]
	Gamma = 0.3333333333;					//	[-]
	rockrho = 2500.0;					//	[kg/m3]
	gravity = 9.81;					//	[m/s2]
	visc =0.001;
	Po = SigTop*Gamma;
	
	segtime = 0.0;						//	[s]
	
	
	segtime = timestep;		
	if (segtime == 0.0) {
		segtime = 1.0e-10;
	}
	
//	int in;
	REAL pD = 0.0, uDx = 0.0,uDy = 0.0, sigDy = 0.0, qDy = 0.0;
	
	sol[0]=0.0;
	sol[1]=0.0;
	sol[2]=0.0;
	
	REAL H=1.0;
	REAL tD = segtime;
	REAL yD = abs(y)/H;
	REAL xi = (yD/(2.0*sqrt(tD)));
	
#ifdef USING_BOOST
	pD = boost::math::erf(xi);
	uDy = -(((2.0*sqrt(tD)*exp(-pow(xi,2)))/(sqrt(PI))) - yD*(boost::math::erfc(xi)));
#else
	pD = erf(xi);
#ifndef WIN32
	uDy = -(((2.0*sqrt(tD)*exp(-pow(xi,2)))/(sqrt(PI))) - yD*(erfc(xi)));
#else
	std::cout << "It is necessary get a implementation of the erfc function." << std::endl;
#endif
#endif
	qDy = -((exp(-pow(xi,2)))/(sqrt(tD)*sqrt(PI)));

	sol[0] = pD;
	sol[2] = uDx;	
	sol[2] = uDy;
	flux(0,0) = 0.0;
	flux(1,0) = qDy;	
	flux(2,0) = 0.0;	
	
	//		Secondary Variables		
	flux(0,0)=0.0; // SigX
	flux(1,0)=sigDy; // SigY
	flux(2,0)=0.0; // TauXY
	flux(3,0)=0.0; // Qx
	flux(4,0)=qDy; // Qy			
	
}


	void ExactSolution2DLineSource(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
		
		// Variable definition
		REAL x = ptx[0];							//	[m]	
		REAL y = ptx[1];							//	[m]	
		REAL z = 0.0;								//	[m]	
		REAL r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));	//	[m]		
		REAL Currenttime = 0.0;							//	[s]		
		if ( r == 0.0) {
			x=1.0e-8;
			y=1.0e-8;
			r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
		}
		
		Currenttime = timestep;
		if (Currenttime == 0.0) {
			Currenttime = 1.0e-8;
		}		
		
		// Parameters definition from: Altmann-J.B._Poroelastic-contribution-to-the-reservoir-stress-path_2010
        REAL lamb, lambu, alpha, G, rhof, S, Gamma;
		lamb = 8.3e9;						//	[Pa]
		lambu = 13.4e9;					//	[Pa]
		alpha = 0.7;						//	[-]
		G= 5.5e9;							//	[Pa]
		rhof = 1000.0;						//	[kg/m3]
		S = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[1/Pa]
		Gamma = ((lambu-lamb)/(alpha*(lambu+2.0*G)));

		REAL qRate, c, kovermu;
        qRate = -20.0;						//	[kg/s]
		c= 0.083;							//	[m2/s]
		kovermu = c*S;
		REAL PI = atan(1.)*4.;

		sol[0]=0.;		// Pressure
		sol[1]=0.;		// Ux
		sol[2]=0.;		// Uy
		flux(0,0)=0.0; // SigX
		flux(1,0)=0.0; // SigY
		flux(2,0)=0.0; // TauXY
		flux(3,0)=0.0; // Qx
		flux(4,0)=0.0; // Qy		
		
		
		REAL Pressure = 0.0;					//	[Pa]
		REAL Ux = 0.0;							//	[m]
		REAL Uy = 0.0;							//	[m]
		REAL Sigx = 0.0;						//	[Pa]
		REAL Sigy = 0.0;						//	[Pa]
		REAL Tauxy = 0.0;						//	[Pa]
		REAL Qx = 0.0;							//	[kg/s]
		REAL Qy = 0.0;							//	[kg/s]		
		
//		Use this Block for normal calculations
//		REAL Zz		= (pow(r, 2)/(4.0*c*Currenttime));
//		Pressure	= (qRate/(4*PI*rhof*kovermu))*Exponential_Integral_Ei(-Zz);
//		Ux			= ((-qRate*alpha*x)/(8*PI*rhof*kovermu*(lamb+2.0*G)))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
//		Uy			= ((-qRate*alpha*y)/(8*PI*rhof*kovermu*(lamb+2.0*G)))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
//		Sigx		= (qRate*alpha*G/(4*PI*rhof*kovermu*(2.0*G+lamb)))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(x,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
//		Sigy		= (qRate*alpha*G/(4*PI*rhof*kovermu*(2.0*G+lamb)))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(y,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
//		Tauxy		= (2.0*qRate*alpha*G*x*y/(4*PI*rhof*kovermu*(2.0*G+lamb)*pow(r,2)))*((1/Zz)*(1-exp(-Zz)));

//		Use this Block for dimensionless calculations		
		REAL Zz = (pow(r,2)/(4*Currenttime));				
		Pressure = -(1/(4*PI))*Exponential_Integral_Ei(-Zz);
		Ux = x*(1/(8*PI))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
		Uy = y*(1/(8*PI))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
		Sigx = (1/(4*PI))*((alpha*G)/(lamb+2.0*G))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(x,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
		Sigy = (1/(4*PI))*((alpha*G)/(lamb+2.0*G))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(y,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
		Tauxy = (1/(4*PI*pow(r,2)))*((2*alpha*G)/(lamb+2.0*G))*(x*y)*((1/Zz)*(1-exp(-Zz)));	
		Qx = (x/(2.0*PI*pow(r,2)))*exp(-Zz);
		Qy = (y/(2.0*PI*pow(r,2)))*exp(-Zz);		

		
//		State Variables
		sol[0] = Pressure;
		sol[1] = Ux*Gamma;
		sol[2] = Uy*Gamma;
		
//		Secondary Variables		
		flux(0,0)=Sigx; // SigX
		flux(1,0)=Sigy; // SigY
		flux(2,0)=Tauxy; // TauXY
		flux(3,0)=Qx; // Qx
		flux(4,0)=Qy; // Qy			
	}

	void SolutionExactRosa1D(TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux) {
		//REAL x = ptx[0];
		REAL x = ptx[0];//-12500;

		//	Parameters
		REAL Eyoung, poisson, alpha;
        Eyoung = 1.43e10;					//	[Pa]
		poisson = 0.3;						//	[-]
		alpha=0.0;///((1.19667e10)*3);//1.70667e10		1.19667e10					//	[-]
		
		REAL Phi= 1.0;//9.60784e-11;//1.35695e-10				//	[-]
        REAL Ct, Se,Visc, Kmed, Qo, Bo, segtime;
		Ct = 1.0;					//	[kg/m3]
		Se = Ct*Phi;							//	[m2/s]
		Visc = 1.0;						//	[Pa.s]
		Kmed = 1.0;//7.8784288e-15;//1.109542e-14						//	[m2]
		Qo = 2.0;
		Bo = 1.0;
		segtime = 0.0;					//	[s]
		
		REAL PI = atan(1.)*4.;
		
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
		REAL Pressure = 0.0;
		#ifdef USING_BOOST
			Pressure = ((0.5*Qo*Bo*Visc)/(Kmed*1.0))*(sqrt((4*Eta*segtime)/PI)*(exp(-(pow(x,2)/(4*Eta*segtime))))-(x*boost::math::erfc(x/sqrt(4*Eta*segtime))));
		#else
#ifndef WIN32
			Pressure = ((0.5*Qo*Bo*Visc)/(Kmed*1.0))*(sqrt((4*Eta*segtime)/PI)*(exp(-(pow(x,2)/(4*Eta*segtime))))-(x*erfc(x/sqrt(4*Eta*segtime))));
#else
	std::cout << "It is necessary get a implementation of the erfc function." << std::endl;
#endif
		#endif
		
		sol[0] = Pressure;
		//	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
		//	sol[2] = (-1.+ sigD)*pini;
	}

	// Analytical solution for flamant problem Elasticity Theory: Applications and Numerics 
	void ExactSolutionFlamantProblem(TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux)
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


	// Exact Solution Madels problem ref	
	void ExactSolutionMandelsProblemwitheffect(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux)
	{
		// Defition of variables
        REAL x, y, z, PI;
		x = ptx[0];
		y = ptx[1];
		z = 0.0;
		PI = atan(1.)*4.;
        REAL B, vu, v, F, L, c;
		B = 0.8;
		vu = 0.4;
		v = 0.2;
		F = 100000.0;
		L = 100.0;
		c = 0.000777778;
		
		sol[0]=0.;
		sol[1]=0.;
		sol[2]=0.;		
		
//		REAL lambdai[14] = {1.28734, 4.6316, 7.80598, 10.9614, 14.1106, 20.402, 23.546, 29.8326,36.1179, 39.2604, 51.829, 58.113, 67.5387, 212.056};
//		REAL lambdai[12] =	{1.45689, 4.67677, 7.83271, 10.9804, 14.1254, 20.4122, 26.6973,51.8331, 67.5418, 70.6835, 89.5335, 105.242};
		
//		REAL lambdai[55] =		{1.45689, 4.67677, 7.83271, 10.9804, 14.1254, 20.4122, 26.6973,
//			29.8395, 39.2657, 45.5494, 51.8331, 61.2583, 64.4001, 67.5418, 
//			70.6835, 76.9669, 83.2502, 89.5335, 105.242, 111.525, 120.95, 
//			130.375, 133.516, 136.658, 161.791, 168.074, 180.641, 190.065, 
//			199.49, 202.632, 218.34, 227.765, 230.906, 234.048, 240.331, 243.473, 
//			256.039, 262.322, 268.606, 274.889, 328.296, 334.579, 350.287, 
//			353.429, 369.137, 372.278, 381.703, 416.261, 491.659, 551.349, 
//			563.916, 664.447, 1091.7, 1141.97, 2863.56};
		
//		REAL lambdai[22] = 	{1.39325, 4.65878, 7.82203, 10.9728, 14.1195, 17.2643, 20.4081,
//			23.5513, 26.6942, 32.9791, 39.2635, 42.4056, 45.5476, 61.257,
//			67.5405, 86.3909, 92.6743, 139.799, 149.224, 183.782, 234.048,
//			413.119};

		REAL lambdai[639] = {1.39325, 4.65878, 7.82203, 10.9728, 14.1195, 17.2643, 20.4081, \
			23.5513, 26.6942, 32.9791, 36.1214, 39.2635, 42.4056, 45.5476, \
			48.6896, 58.1152, 61.257, 64.3988, 67.5405, 70.6823, 73.824, 80.1075, \
			83.2492, 86.3909, 89.5326, 92.6743, 102.099, 105.241, 111.524, \
			114.666, 120.949, 124.091, 127.233, 130.374, 133.516, 136.657, \
			139.799, 146.082, 149.224, 152.366, 155.507, 164.932, 168.074, \
			171.215, 174.357, 177.499, 183.782, 190.065, 193.207, 196.348, \
			202.631, 205.773, 212.056, 215.198, 218.34, 227.764, 234.048, \
			237.189, 240.331, 246.614, 252.897, 256.039, 259.18, 262.322, \
			265.464, 268.605, 271.747, 274.888, 278.03, 281.172, 284.313, \
			287.455, 290.596, 293.738, 296.88, 300.021, 303.163, 306.304, \
			309.446, 315.729, 318.871, 322.012, 325.154, 328.296, 331.437, \
			334.579, 340.862, 344.004, 347.145, 350.287, 353.428, 362.853, \
			365.995, 369.136, 372.278, 375.42, 378.561, 381.703, 384.844, \
			391.128, 394.269, 397.411, 403.694, 406.836, 409.977, 413.119, \
			416.26, 422.544, 425.685, 428.827, 435.11, 438.252, 441.393, 444.535, \
			450.818, 453.96, 457.101, 463.384, 466.526, 469.668, 472.809, \
			479.092, 482.234, 485.376, 491.659, 501.084, 504.225, 507.367, \
			513.65, 519.933, 523.075, 526.216, 529.358, 532.499, 535.641, \
			538.783, 541.924, 545.066, 548.207, 551.349, 554.491, 557.632, \
			560.774, 563.915, 567.057, 570.199, 576.482, 582.765, 585.907, \
			589.048, 592.19, 595.331, 598.473, 601.615, 604.756, 607.898, \
			611.039, 614.181, 620.464, 623.606, 626.747, 629.889, 633.031, \
			636.172, 639.314, 645.597, 648.738, 651.88, 655.022, 658.163, \
			661.305, 664.446, 667.588, 670.73, 673.871, 677.013, 680.154, \
			683.296, 686.438, 689.579, 692.721, 695.862, 699.004, 702.146, \
			705.287, 708.429, 714.712, 720.995, 724.137, 727.278, 730.42, \
			733.562, 736.703, 739.845, 742.986, 746.128, 749.27, 752.411, \
			755.553, 758.694, 764.977, 768.119, 771.261, 774.402, 777.544, \
			780.685, 783.827, 786.969, 790.11, 793.252, 796.393, 799.535, \
			802.677, 808.96, 812.101, 815.243, 818.385, 821.526, 824.668, \
			827.809, 830.951, 834.093, 837.234, 840.376, 852.942, 856.084, \
			859.225, 862.367, 865.508, 868.65, 871.792, 878.075, 881.216, \
			884.358, 890.641, 896.924, 900.066, 903.208, 906.349, 912.632, \
			915.774, 918.916, 922.057, 925.199, 928.34, 931.482, 934.624, \
			944.048, 950.332, 953.473, 956.615, 959.756, 962.898, 966.039, \
			969.181, 972.323, 978.606, 984.889, 988.031, 991.172, 994.314, \
			997.455, 1000.6, 1006.88, 1010.02, 1013.16, 1025.73, 1032.01, \
			1044.58, 1047.72, 1050.86, 1054., 1057.15, 1063.43, 1069.71, 1072.85, \
			1082.28, 1085.42, 1088.56, 1091.7, 1094.84, 1097.99, 1110.55, \
			1119.98, 1123.12, 1132.54, 1135.69, 1138.83, 1141.97, 1145.11, \
			1148.25, 1151.39, 1154.54, 1163.96, 1173.38, 1176.53, 1179.67, \
			1182.81, 1192.23, 1201.66, 1204.8, 1207.94, 1211.08, 1223.65, \
			1226.79, 1229.93, 1233.07, 1236.22, 1239.36, 1242.5, 1245.64, \
			1248.78, 1264.49, 1267.63, 1270.77, 1273.92, 1277.06, 1280.2, \
			1286.48, 1308.47, 1321.04, 1324.18, 1330.46, 1333.61, 1336.75, \
			1339.89, 1343.03, 1349.31, 1368.16, 1374.45, 1390.15, 1393.3, \
			1399.58, 1409., 1412.15, 1415.29, 1427.85, 1434.14, 1437.28, 1440.42, \
			1446.7, 1452.99, 1462.41, 1468.69, 1471.84, 1478.12, 1496.97, \
			1506.39, 1515.82, 1522.1, 1537.81, 1547.23, 1550.38, 1559.8, 1562.94, \
			1566.08, 1572.37, 1578.65, 1581.79, 1588.07, 1594.36, 1597.5, \
			1613.21, 1638.34, 1641.48, 1657.19, 1660.33, 1663.47, 1666.61, \
			1679.18, 1691.75, 1694.89, 1710.6, 1726.31, 1729.45, 1732.59, \
			1735.73, 1742.01, 1757.72, 1760.86, 1767.15, 1773.43, 1792.28, \
			1795.42, 1807.99, 1817.41, 1826.84, 1836.26, 1839.4, 1845.69, \
			1861.39, 1867.68, 1873.96, 1877.1, 1889.67, 1892.81, 1902.23, \
			1905.38, 1930.51, 1936.79, 1946.22, 1965.07, 1990.2, 1993.34, \
			1996.48, 2002.77, 2015.33, 2018.47, 2024.76, 2037.32, 2049.89, \
			2062.46, 2068.74, 2071.88, 2084.45, 2097.01, 2115.86, 2119., 2144.14, \
			2162.99, 2175.55, 2191.26, 2206.97, 2213.25, 2216.39, 2238.38, \
			2247.81, 2272.94, 2285.51, 2288.65, 2301.22, 2316.92, 2329.49, \
			2335.77, 2345.2, 2348.34, 2373.47, 2382.9, 2389.18, 2392.32, 2395.46, \
			2398.61, 2411.17, 2430.02, 2448.87, 2452.01, 2467.72, 2483.43, \
			2508.56, 2533.69, 2543.12, 2552.54, 2555.69, 2561.97, 2565.11, \
			2574.54, 2590.24, 2596.53, 2602.81, 2671.92, 2678.21, 2684.49, \
			2700.2, 2703.34, 2706.48, 2709.62, 2712.77, 2725.33, 2737.9, 2753.61, \
			2775.6, 2778.74, 2825.86, 2829., 2835.29, 2838.43, 2841.57, 2857.28, \
			2879.27, 2891.84, 2901.26, 2907.54, 2929.54, 2932.68, 2948.38, \
			2960.95, 2970.38, 2992.37, 2995.51, 3004.93, 3008.07, 3058.34, \
			3064.62, 3070.91, 3092.9, 3105.46, 3111.75, 3124.31, 3140.02, \
			3155.73, 3168.3, 3209.14, 3218.56, 3237.41, 3243.69, 3256.26, \
			3275.11, 3293.96, 3319.09, 3331.66, 3337.94, 3363.07, 3369.36, \
			3385.07, 3397.63, 3422.77, 3435.33, 3460.46, 3491.88, 3501.3, \
			3507.59, 3513.87, 3523.3, 3526.44, 3529.58, 3539., 3561., 3592.41, \
			3608.12, 3623.83, 3626.97, 3661.53, 3667.81, 3692.94, 3714.93, \
			3733.78, 3736.92, 3799.76, 3821.75, 3843.74, 3856.3, 3881.44, \
			3912.85, 3922.28, 4010.24, 4035.38, 4060.51, 4063.65, 4069.93, \
			4076.22, 4079.36, 4117.06, 4126.48, 4151.61, 4205.02, 4208.16, \
			4230.15, 4261.57, 4274.14, 4286.7, 4327.54, 4358.96, 4362.1, 4365.24, \
			4387.23, 4459.49, 4490.91, 4497.19, 4541.17, 4556.88, 4685.69, \
			4720.24, 4723.38, 4735.95, 4792.5, 4839.62, 4842.77, 4855.33, \
			4874.18, 4883.61, 4905.6, 4915.02, 5043.83, 5075.24, 5112.94, \
			5213.47, 5298.3, 5357.99, 5367.41, 5442.81, 5502.5, 5537.06, 5631.3, \
			5697.28, 5709.84, 5782.1, 5882.63, 5967.46, 5998.87, 6033.43, \
			6049.14, 6121.39, 6124.53, 6133.96, 6137.1, 6203.07, 6303.61, 6454.4, \
			6602.06, 6646.04, 6884.8, 6922.5, 7230.38, 7330.91, 7415.73, 7937.23, \
			7990.64, 8116.3, 8125.73, 8150.86, 8179.14, 8320.51, 8383.34, \
			8424.18, 8527.85, 8568.69, 8835.73, 9401.22, 9520.6, 9552.01, \
			9668.25, 9690.24, 10061., 10516.5, 10582.5, 10767.8, 11229.6, \
			13997.4, 17283.5};	

		REAL H, pD, qD, uxD, uyD, tD;
        H=1.0;
		pD = 0.0;
		qD = 0.0;
		uxD = 0.0;
		uyD = 0.0;
		tD = timestep;
//		REAL tD = timestep*(c/pow(L,2));		
//		REAL xD = x/L;
		REAL xD = x;
		REAL yD = y;		
		for (int in =0; in < 639; in++) 
		{

			pD += ((sin(lambdai[in]))/(lambdai[in]-sin(lambdai[in])*cos(lambdai[in])))*(cos(lambdai[in]*xD)-cos(lambdai[in]))*(exp(-pow(lambdai[in],2)*tD));
			qD += ((lambdai[in]*sin(lambdai[in]))/(lambdai[in]-sin(lambdai[in])*cos(lambdai[in])))*(sin(lambdai[in]*xD))*(exp(-pow(lambdai[in],2)*tD));
			uxD += xD*((v-vu)/2)*((sin(lambdai[in]))/(lambdai[in]-sin(lambdai[in])*cos(lambdai[in])))*(cos(lambdai[in]))*(exp(-pow(lambdai[in],2)*tD))+
				((cos(lambdai[in]))/(lambdai[in]-sin(lambdai[in])*cos(lambdai[in])))*(sin(lambdai[in]*xD))*(exp(-pow(lambdai[in],2)*tD));
			uyD += yD*(-((1-v)/2)+((1-vu)))*((sin(lambdai[in]))/(lambdai[in]-sin(lambdai[in])*cos(lambdai[in])))*(cos(lambdai[in]))*(exp(-pow(lambdai[in],2)*tD));
					
//			sigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
		}		
		
//		sol[0] = pD*((2.*F*B*(1+vu))/(3*L));
		sol[0] = pD*2.0;
//		sol[1] = uxD;
//		sol[2] = uyD;		
		
		flux(0,0) = 0.0;
		flux(2,0) = 0.0;			
//		flux(1,0) = qD*((2.0*B*(1+vu))/(3));		
		
	}

	void SolutionExactRosa1DPseudo(TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
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

	void InitialPressureDistribution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol)
	{
		REAL x, y, z;
        x = ptx[0];
		y = ptx[1];
		z = ptx[2];
		sol[0]=0.;
//		sol[1]=0.;
//		sol[2]=0.;
		REAL Pressure = 0.0;		
		sol[0] = Pressure;
		
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
