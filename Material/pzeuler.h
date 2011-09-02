/**
 * \file
 * @brief Contains the TPZEulerEquation class which implements the weak statement of the compressible euler equations.
 */
//$Id: pzeuler.h,v 1.7 2009-12-15 17:28:28 phil Exp $

#ifndef PZEULER_H
#define PZEULER_H

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzdiscgal.h"

#include "pzausmflux.h"
#include "pzgradientflux.h"
#include "pzlog.h"
class TPZCompMesh;

/**
 * @ingroup material
 * @brief This material implements the weak statement of the compressible euler equations \n
 * for Olivier Roussel's project.
 */
/**
 * It is for transient analysis, finite volume method and
 * explicit time integrator only.
 */
class TPZEulerEquation : public TPZDiscontinuousGalerkin{
	
public:
	
	enum BCType{EFreeSlip = 1};
	
	enum CALCType{EFlux = 1, EGradient = 2};
	
	static void SetComputeFlux(){
		gType = EFlux;
	}
	
	static void SetComputeGradient(){
		gType = EGradient;
	}
	
#ifdef LinearConvection
	static void SetLinearConvection(TPZCompMesh * cmesh, TPZVec<REAL> &Celerity);
#endif
	
private:
	
#ifdef LinearConvection
	TPZVec<REAL> fCelerity;
#endif
	
	static CALCType gType;
	
	/** @brief Ratio between specific heat is constant and the specific heat the constant volume of a polytropic gas */
	static REAL gGamma;
	
	/** @brief Convective flux object */
	TPZAUSMFlux fAUSMFlux;
	
	/** @brief Gradient flux object */
	TPZGradientFlux fGradientFlux;
	
	/** @brief Compute Euler Flux */
	void ComputeEulerFlux(TPZVec<REAL> &sol, TPZFMatrix & F);
	
public:
	
	static REAL Gamma(){ return gGamma; }
	
	/** @brief Convert from primitive to conservative variables */
	static void FromPrimitiveToConservative(TPZVec<REAL> &sol,REAL gamma);
	
	/** @brief Convert from conservative to primitive variables */
	static void FromConservativeToPrimitive(TPZVec<REAL> &sol,REAL gamma);
	
	/** @brief Constructor with Gamma value */
	TPZEulerEquation(int nummat, REAL gamma);
	
	/** @brief Default destructor */
	~TPZEulerEquation();
	
	/** @brief Default constructor */
	TPZEulerEquation();
	
	/** @brief Copy constructor */
	TPZEulerEquation(const TPZEulerEquation &cp);
	
	/** @brief Creates a copy of this */
	TPZAutoPointer<TPZMaterial> NewMaterial();
	
	/** @brief Object-based overload */
	int NStateVariables();
	
	/** @brief Object-based overload */
	virtual int Dimension();
	
	/** @brief Returns the pressure value */
	static REAL Pressure(TPZVec<REAL> &U, double gamma);
	
	/** @brief Computes sound speed */
	REAL cSpeed(TPZVec<REAL> & sol);
	
	/** @brief Returns \f$ u = Sqrt(u2 + v2 + w2) \f$ */
	REAL uRes(TPZVec<REAL> & sol);
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name(){return "TPZEulerEquation";}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

	/** 
	 * @name Contribute methods 
	 * @brief data contains material data. data.soll and data.solr are expected in primitive variables
	 */
	/** @{ */
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek, TPZFMatrix &ef,
							  TPZBndCond &bc);
	
	virtual void ContributeBCInterface(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ek,TPZFMatrix &ef,
									   TPZBndCond &bc);
	
	virtual void ContributeBCInterface(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc);
	/** @} */
	
};

#endif///PZEULER_H

