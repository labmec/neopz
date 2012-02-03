/*
 *  Elasticitytensordata.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 2/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "Elasticitytensordata.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
	
    /**
     * Void Constructor
     */	
	ElasticityTensorData::ElasticityTensorData()
	{
		//	-999 means no initialized
		
		fCelastic.Resize(6,6);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				fCelastic(i,j)=-999.0;
			}
		}		
	fvp = -999.0;
	fvs = -999.0;
	fmu = -999.0;
	flambda = -999.0;
	fEyoung = -999.0;
	fPnu = -999.0;	
		
				
	}
	

	/**
	 * Destrutor
	 */
	ElasticityTensorData::~ElasticityTensorData()
	{
	}
	
	//  Set Dimension and kind of problem 
	//	1D, 2D, 3D
	//	Elastic linear strain 1D, kindflag = 0
	//	Elastic linear stress 1D, kindflag = 1
	//	Elastic plane  strain 2D, kindflag = 0
	//	Elastic plane  strain 2D, kindflag = 1
	//	Elastic 3D	
	
	
// Not check this conversions, recommendation: use SetDataLame
	
	//  Set Data for elasticity tensor relation from rock P and S wave velocities vp, vs and rock density 
	void ElasticityTensorData::SetDataSeismic(double vp, double vs, double rockrho)
	{
		fvp = vp;
		fvs = vs;
		fmu = rockrho*vs*vs;
		flambda = rockrho*((vp*vp)-2*(vs*vs)); // note:use this with Vp > 2*Vs
	}
	
	//  Set Data for elasticity tensor relation from Young and Poisson modulus 
	void ElasticityTensorData::SetDataHook(double Eyoung, double Pnu)
	{
		fEyoung = Eyoung;
		fPnu = Pnu;
		fmu = (Eyoung)/(2*(1+Pnu));
		flambda = (Eyoung*Pnu)/((1+Pnu)*(1-2*Pnu));
	}
	
	//  Set Data for elasticity tensor relation from lamé parameters 
	void ElasticityTensorData::SetDataLame(double mu, double lambda)
	{
		fmu = mu;
		flambda = lambda;
	}
	
	//	Set kind of problem
	// Theory for this implementation can found on pag. 69 Linear Elastic models
	// link: http://geodynamics.org/cig/software/pylith/pylith_manual-1.6.1.pdf
	void ElasticityTensorData::Setproblem(double Dimension, double kindflag)
	{
		fDimension = Dimension;  
		fkindflag = kindflag;
		
		if (Dimension == 1 && kindflag == 0) {
			// Linear strain
			fCelastic.Resize(1, 1);
			fCelastic(0,0) = (flambda + 2*fmu);
		}
			else {
				if (Dimension == 1 && kindflag == 1) {
					// Linear stress
					fCelastic.Resize(1, 1);
					fCelastic(0,0) = (fmu*(3*flambda + 2*fmu))/(flambda + fmu);
				}
					else {
						if (Dimension == 2 && kindflag == 0) {
							// Plane strain
							fCelastic.Resize(3, 3);
							fCelastic(0,0) = flambda + 2*fmu; 
							fCelastic(1,0) = flambda;
							fCelastic(2,0) = 0.0;
							fCelastic(0,1) = flambda;
							fCelastic(1,1) = flambda + 2*fmu;
							fCelastic(2,1) = 0.0;
							fCelastic(0,2) = 0.0;
							fCelastic(1,2) = 0.0;
							fCelastic(2,2) = 2*fmu;
						}
							else {
								if (Dimension == 2 && kindflag == 1) {
									// Plane Stress
									fCelastic.Resize(3, 3);
									fCelastic(0,0) = (4*fmu*flambda)/(flambda + 2*fmu); 
									fCelastic(1,0) = (4*fmu*(flambda + fmu))/(flambda + 2*fmu);
									fCelastic(2,0) = 0.0;
									fCelastic(0,1) = (4*fmu*flambda)/(flambda + 2*fmu);
									fCelastic(1,1) = (4*fmu*(flambda + fmu))/(flambda + 2*fmu);
									fCelastic(2,1) = 0.0;
									fCelastic(0,2) = 0.0;
									fCelastic(1,2) = 0.0;
									fCelastic(2,2) = 2*fmu;
								}
									else {
										if (Dimension == 3) {
											// General 3D
											fCelastic.Resize(6, 6);
											fCelastic(0,0) = flambda + 2*fmu; // C1111
											fCelastic(1,0) = flambda; // C1122
											fCelastic(2,0) = flambda; // C1133
											fCelastic(3,0) = 0.0; // C1112
											fCelastic(4,0) = 0.0; // C1123
											fCelastic(5,0) = 0.0; // C1113
											fCelastic(0,1) = flambda; // C1122
											fCelastic(1,1) = flambda + 2*fmu; // C2222
											fCelastic(2,1) = flambda; // C2233
											fCelastic(3,1) = 0.0; // C2212
											fCelastic(4,1) = 0.0; // C2223
											fCelastic(5,1) = 0.0; // C2213
											fCelastic(0,2) = flambda; // C1133
											fCelastic(1,2) = flambda; // C2233
											fCelastic(2,2) = flambda + 2*fmu; // C3333
											fCelastic(3,2) = 0.0; // C3312
											fCelastic(4,2) = 0.0; // C3323
											fCelastic(5,2) = 0.0; // C3313
											fCelastic(0,3) = 0.0; // C1112
											fCelastic(1,3) = 0.0; // C2212
											fCelastic(2,3) = 0.0; // C3312
											fCelastic(3,3) = 2*fmu; // C1212
											fCelastic(4,3) = 0.0; // C1223
											fCelastic(5,3) = 0.0; // C1213
											fCelastic(0,4) = 0.0; // C1123
											fCelastic(1,4) = 0.0; // C2223
											fCelastic(2,4) = 0.0; // C3323
											fCelastic(3,4) = 0.0; // C1223
											fCelastic(4,4) = 2*fmu; // C2323
											fCelastic(5,4) = 0.0; // C2313
											fCelastic(0,5) = 0.0; // C1113
											fCelastic(1,5) = 0.0; // C2213
											fCelastic(2,5) = 0.0; // C3313
											fCelastic(3,5) = 0.0; // C1213
											fCelastic(4,5) = 0.0; // C2313
											fCelastic(5,5) = 2*fmu; // C1313
										}
											else {
												PZError << "\n inconsistent input data (Dimension 1,2,3 or kindflag 0,1): \n" << "Dimension = " << fDimension << " fkinflag = " << fkindflag << "\n";
												return;
											}
									}
							}
					}
			}
	}

	

