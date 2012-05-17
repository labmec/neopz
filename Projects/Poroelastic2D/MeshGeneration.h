/*
 *  MeshGeneration.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/15/12.
 *  Copyright 2012 __LabMeC__. All rights reserved.
 *
 */

//---------------------------------------------------------------------------

#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzgmesh.h"

#ifndef MeshGenerationH
#define MeshGenerationH
//---------------------------------------------------------------------------


class MeshGeneration
{
	
public:
	
    /**
     * Void Constructor
     */	
	MeshGeneration();
	
    /**
     * Destrutor
     */
	~MeshGeneration();
	
	//	Set kind of problem
	void Setdimensions();	
	
	void Setdimensions(REAL xdimension, REAL ydimension, REAL zThickness);

	void Setdimensions(REAL rdimension, REAL zThickness);
	
	
	//  2D Geometric Mesh Generation 
	
	// Geometries for differents Problems
	
	//	Set kind of Geometry for 1D Validation
	/** @brief Cretaes Rectangular Grid for 1D Validation Problem
	 * 1D Validation Problem is given on: */	
	TPZGeoMesh * GeometricMesh1DValidation(double Dimension, double kingflag);	
	
	TPZGeoMesh	* MalhaGeom(TPZVec < int > matIdlist);	
	
	//	Set kind of Geometry for 2D Validation
	/** @brief Cretaes Circular Grid for 2D Validation Problem
	 * 2D Validation Problem is given on: Rudnicki 1986 */
	TPZGeoMesh * GeometricMesh2DValidation(TPZVec <int> matIdlist);
	
	
	TPZGeoMesh * MalhaGeoGravenobj(int nLayers, REAL LlengthFootFault, REAL DipFaultAngleleft, REAL DipFaultAngleright, REAL WellFaultlength, TPZVec <bool> Productionlayer, bool InterfaceElement);
	
//	//  Set Data for elasticity tensor relation from Young and Poisson modulus 
//	void GeometricMesh1DValidation(double Eyoung, double Pnu);
//	
//	//  Set Data for elasticity tensor relation from lamé parameters 
//	void GeometricMesh1DValidation(double mu, double lambda);
	
	//  Recovery material properties
//	double getVp() { return fvp; }
//	double getVs() { return fvs; }
//	double getEyoug() { return fEyoung; }
//	double getPnu() { return fPnu; }	
//	double getMu() { return fmu; }
//	double getLambda() { return flambda; }	
//	double getRockRho() {return frockrho; }
	
protected:
	
	double fDimension, fkindflag, fEyoung, fPnu, fmu, flambda, fvp, fvs, frockrho;
	double fxLength, fyLength, frLength, fNnodes, fThickness;	
    /**
     * Constant elastic tensor C
     */
    TPZFMatrix<REAL> fCelastic;
	
	/**
     * Vector mapping m
     */
	TPZFMatrix<REAL> fm;	
};

#endif