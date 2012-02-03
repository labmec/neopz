/*
 *  Elasticitytensordata.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 2/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

//---------------------------------------------------------------------------

#include "pzfmatrix.h"
#include "pzmatrix.h"

#ifndef ElasticityTensorDataH
#define ElasticityTensorDataH
//---------------------------------------------------------------------------


class ElasticityTensorData
{
	
public:
	
    /**
     * Void Constructor
     */	
	ElasticityTensorData();
	
    /**
     * Destrutor
     */
	~ElasticityTensorData();

//  Set Dimension and kind of problem 
//	1D, 2D, 3D
//	Elastic linear strain 1D, kindflag = 0
//	Elastic linear stress 1D, kindflag = 1
//	Elastic plane  strain 2D, kindflag = 0
//	Elastic plane  strain 2D, kindflag = 1
//	Elastic 3D	

	
//	Set kind of problem
	void Setproblem(double Dimension, double kingflag);	
	
//  Set Data for elasticity tensor relation from rock P and S wave velocities vp, vs and rock density 
	void SetDataSeismic(double vp, double vs, double rockrho);

//  Set Data for elasticity tensor relation from Young and Poisson modulus 
	void SetDataHook(double Eyoung, double Pnu);
	
//  Set Data for elasticity tensor relation from lamé parameters 
	void SetDataLame(double mu, double lambda);

//  Recovery material properties
	double getVp() { return fvp; }
	double getVs() { return fvs; }
	double getEyoug() { return fEyoung; }
	double getPnu() { return fPnu; }	
	double getMu() { return fmu; }
	double getLambda() { return flambda; }	
	double getRockRho() {return frockrho; }
	
protected:
	
//Not Complete the documentation yet
	
	/**
	 * fvp: P Wave velocity  [m/s]
	 * fvs: S Wave velocity  [m/s]
	 * frockrho: Rock density [kg/m3]
	 * fEyoung: Young's modulus -> units may be [Pa]
	 * fPnu: Poisson modulus	 
	 * fmu: Lame parameter
	 * flambda: Lame parameter
	 * fDimension: Problem dimenion 
	 * fkindflag; strain (0) or stress (1) problem 
	 */
	double fDimension, fkindflag, fEyoung, fPnu, fmu, flambda, fvp, fvs, frockrho;
    /**
     * Constant elastic tensor C
     */
    TPZFMatrix fCelastic;

	/**
     * Vector mapping m
     */
	TPZFMatrix fm;	
};

#endif
