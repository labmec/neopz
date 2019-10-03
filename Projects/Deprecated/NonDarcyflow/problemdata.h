#ifndef TPZProblemDATAH
#define TPZProblemDATAH
/*
 *  problemdata.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include "pzdiscgal.h"
#include <math.h>

#ifdef _AUTODIFF
#include "fad.h"

/** @brief State: Stiffness or Mass Matrix Calculations */
enum EState { ELastState = 0, ECurrentState = 1 };

typedef TFad<3, REAL> VarFad;

class TPZProblemDATA {
	
public:
	

	EState fState;
	
	REAL fRhoRef;
	
	REAL fViscRef;
	
	REAL fPorRef;

	REAL fCRho;
	
	REAL fCVisc;
	
	REAL fCPor;	
	
	REAL fKh;
	
	REAL fKv;
	
	REAL fDeltaT;
	
	REAL fPini;
	
	REAL fNSteps;
	
	REAL fRw;
	
	REAL fRe;
	
	REAL fH;
	
	
	/** @brief Characteristic length - m */
	REAL fLref;	
	
	/** @brief Permeability reference - m2 */
	REAL fKref;
	
	/** @brief Pressure reference - Pa */	
	REAL fPref;
	
	/** @brief density reference - kg/m3 */
	REAL fRhoref;
	
	/** @brief viscosity reference - Pa s */
	REAL fEtaref;
	
	
	/** @brief K map */	
	TPZStack< TPZFMatrix<REAL> > fKabsoluteMap;	
	
	/** @brief absolute permeability */	
	TPZFMatrix<REAL> fKab;
	
	/** @brief Use or not K map */		
	bool fYorN;
	
	/** @brief Parameter representing temporal scheme for conservation equation */
	REAL fGamma;
	
public:
	
	TPZProblemDATA();
	
	~TPZProblemDATA();
//
//	
//	/** 
//	 * @brief \f$ Rock porosity. \f$ Phi = Phi( p ) \f$
//	 * @param po Refrence pressure
//	 */	
//	void Porosity(REAL po, REAL &poros, REAL &dPorosDp);
//	
//	/** 
//	 * @brief \f$ Oil density RhoOil = RhoOil( po ) \f$
//	 * @param po Refrence pressure
//	 */
//	void RhoOil(REAL po, REAL &RhoOil, REAL &dRhoOilDpo);
//	
//	/** 
//	 * @brief \f$ Water density RhoWater = RhoWater( pw ) \f$
//	 * @param pw Refrence pressure
//	 */
//	void RhoWater(REAL pw, REAL &RhoWater, REAL &dRhoWaterDpo);
//	
//	/** 
//	 * @brief Oil viscosity. \f$ OilViscosity = ViscOleo( po ) \f$
//	 * @param po Refrence pressure
//	 */
//	void OilViscosity(REAL po, REAL &OilViscosity, REAL &dOilViscosityDpo);
//	
//	/** 
//	 * @brief Water viscosity. \f$ WaterViscosity = WaterViscosity( pw ) \f$
//	 * @param po Refrence pressure
//	 */	
//	void WaterViscosity(REAL po, REAL &WaterViscosity, REAL &dWaterViscosityDpo);
//	
//	
//	
//	/** 
//	 * @brief \f$ Rock porosity. \f$ Phi = Phi( p ) \f$
//	 * @param po Refrence pressure
//	 */	
//	void Porosity(VarFad po, VarFad &poros);	
//	
//	/** 
//	 * @brief \f$ Oil density RhoOil = RhoOil( po ) \f$
//	 * @param po Refrence pressure
//	 */
//	void RhoOil(VarFad po, VarFad &RhoOil);	
//	
//	/** 
//	 * @brief \f$ Water density RhoWater = RhoWater( pw ) \f$
//	 * @param pw Refrence pressure
//	 */
//	void RhoWater(VarFad pw, VarFad &RhoWater);		
//	
//	/** 
//	 * @brief Oil viscosity. \f$ OilViscosity = ViscOleo( po ) \f$
//	 * @param po Refrence pressure
//	 */
//	void OilViscosity(VarFad po, VarFad &OilViscosity);	
//	
//	/** 
//	 * @brief Water viscosity. \f$ WaterViscosity = WaterViscosity( pw ) \f$
//	 * @param po Refrence pressure
//	 */	
//	void WaterViscosity(VarFad po, VarFad &WaterViscosity);		
//	
//	
//	
//
//	///Sets permeability data
//	void SetK( TPZFMatrix<STATE> &K ){ this->fKab = K; }
//	
//	/** @brief Set characteristic length - m */
//	REAL SetLreference(REAL &Lref){ fLref = Lref;}
//	
//	/** @brief Set permeability reference - m2 */
//	REAL SetKreference(REAL &Kref){ fKref = Kref;}
//	
//	/** @brief Set pressure reference - Pa */
//	REAL SetPreference(REAL &Pref){ fPref = Pref;}	
//	
//	/** @brief Set density reference - kg/m3 */
//	REAL SetRhoSCreference(REAL &Densityref){ fRhoref = Densityref;}
//	
//	/** @brief Set viscosity reference - Pa s */
//	REAL SetEtaSCreference(REAL &Viscosityref){ fEtaref = Viscosityref;}
//	
//	/** @brief Oil density on standard conditions - kg/m3 */
//	REAL RhoOilSC();
//	
//	/** @brief Water density on standard conditions - kg/m3 */
//	REAL RhoWaterSC();
//	
//	/** @brief Gravity */
//	TPZFMatrix<REAL> Gravity();
//	
//	/** @brief Absolute permeability. */
//	void K(TPZFMatrix<REAL> &Kab);
//	
//	/** @brief Absolute permeability inverse. */
//	TPZFMatrix<REAL>  Kinv(TPZFMatrix<REAL> &Kab);	
//	
//	/** @brief Absolute permeability. */
//	void LoadKMap(std::string MaptoRead);
//	
//	
//	/** @brief Defines simulation time step. */
//	void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}	
//	
//	/** @brief Defines stemporal scheme. */
//	void SetTScheme(REAL timegamma){ this->fGamma = timegamma;}
//	
//	void SetLastState(){ gState = ELastState;}
//	
//	void SetCurrentState(){ gState = ECurrentState;}	
	
	
};

extern TPZProblemDATA globData;



#endif

#endif