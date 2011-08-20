/**
 * \file
 * @brief Contains the TPZAUSMFlux class.
 */
//$Id: pzausmflux.h,v 1.2 2009-08-28 19:43:43 fortiago Exp $

#ifndef TPZAUSMFLUX_H
#define TPZAUSMFLUX_H

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzreal.h"

/**
 * @ingroup material
 * @brief Implements the numerical flux for AUSM problem. (Jorge?)
 */
class TPZAUSMFlux{
	
private:
	
	/**
	 * @brief Ratio between specific heat is constant and the specific heat the constant
	 * volume of a polytropic gas
	 */
	REAL fGamma;
	
	/** @brief Method constants */
	REAL fAlpha, fBeta;
	
public:
	/** @brief Constructor with Gamma value */
	TPZAUSMFlux(REAL gamma);
	/** @brief Copy constructor */
	TPZAUSMFlux(const TPZAUSMFlux &cp);
	
	/** @brief Computes numerical flux */
	void ComputeFlux(TPZVec<REAL> &solL, TPZVec<REAL> &solR, TPZVec<REAL> &normal, TPZVec<REAL> & F);
	
private:
	
	/** @brief Returns sound speed */
	REAL SoundSpeed(TPZVec<REAL> &sol, REAL press);
	
	/** @brief Returns pressure values */
	REAL Pressure(TPZVec<REAL> &sol);
	
	/** @brief Returns speed */
	REAL Speed(TPZVec<REAL> &sol, TPZVec<REAL> &normal, REAL &NormalSpeed);
	
	/** @brief Returns enthalpy */
	REAL Enthalpy(REAL soundSpeed, REAL speed);
	
	/** @brief Auxiliar method only */
	void ComputeInitialData(TPZVec<REAL>&sol,TPZVec<REAL> &normal, REAL&soundSpeed, 
							REAL &Speed, REAL &NormalSpeed, REAL &Enthalpy, REAL &press);
	
	/** @brief Returns pressure in the face */
	REAL FacePressure(REAL pL, REAL pR, REAL Ml, REAL Mr);
	
	/** @brief Returns mach number in the face */
	REAL FaceMachNumber(REAL Ml, REAL Mr);
	
	/** @brief Computes the numerical sound speed at the face */
	REAL NumSoundSpeed(REAL LeftSoundSpeed,REAL RightSoundSpeed);
	
	/** @brief Returns the mass flux */
	REAL MassFlux(REAL NumericalSoundSpeed, REAL rhoL, REAL rhoR, REAL FaceMach);
	
};

#endif
