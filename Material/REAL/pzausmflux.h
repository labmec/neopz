/**
 * \file
 * @brief Contains the TPZAUSMFlux class.
 */

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
	
	/** @brief Ratio between specific heat is constant and the specific heat the constant volume of a polytropic gas */
	STATE fGamma;
	
	/** @brief Method constants */
	STATE fAlpha, fBeta;
	
public:
	/** @brief Constructor with Gamma value */
	TPZAUSMFlux(STATE gamma);
	/** @brief Copy constructor */
	TPZAUSMFlux(const TPZAUSMFlux &cp);
	
	/** @brief Computes numerical flux */
	void ComputeFlux(TPZVec<STATE> &solL, TPZVec<STATE> &solR, TPZVec<REAL> &normal, TPZVec<STATE> & F);
	
private:
	
	/** @brief Returns sound speed */
	STATE SoundSpeed(TPZVec<STATE> &sol, STATE press);
	
	/** @brief Returns pressure values */
	STATE Pressure(TPZVec<STATE> &sol);
	
	/** @brief Returns speed */
	STATE Speed(TPZVec<STATE> &sol, TPZVec<REAL> &normal, STATE &NormalSpeed);
	
	/** @brief Returns enthalpy */
	STATE Enthalpy(STATE soundSpeed, STATE speed);
	
	/** @brief Auxiliar method only */
	void ComputeInitialData(TPZVec<STATE>&sol,TPZVec<REAL> &normal, STATE&soundSpeed, 
							STATE &Speed, STATE &NormalSpeed, STATE &Enthalpy, STATE &press);
	
	/** @brief Returns pressure in the face */
	STATE FacePressure(STATE pL, STATE pR, STATE Ml, STATE Mr);
	
	/** @brief Returns mach number in the face */
	STATE FaceMachNumber(STATE Ml, STATE Mr);
	
	/** @brief Computes the numerical sound speed at the face */
	STATE NumSoundSpeed(STATE LeftSoundSpeed,STATE RightSoundSpeed);
	
	/** @brief Returns the mass flux */
	STATE MassFlux(STATE NumericalSoundSpeed, STATE rhoL, STATE rhoR, STATE FaceMach);
	
};

#endif
