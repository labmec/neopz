//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_pzmultiphase_h
#define PZ_pzmultiphase_h

#include "pzmaterial.h"
#include "pzdiscgal.h"
#ifdef _AUTODIFF
#include "fad.h"
#endif
#include <iostream>
#include <fstream>
#include <string>

/**
 * @ingroup material
 * @author Omar Duran
 * @since 19/08/2013
 * @brief Material to solve a 2d multiphase transport problem by multiphysics simulation
 * @brief Here is used L, Hdiv ... spaces and first order upwind scheme
 */


class TPZMultiphase : public TPZDiscontinuousGalerkin {
    
protected:
	
    /** @brief Problem dimension */
	int fDim;
    
    /** @brief Material id */
    int fmatId;	
	
	/** @brief Definition of constants */
	REAL ff;
	
	/** @brief State: Stiffness or Mass Matrix Calculations */
	enum EState { ELastState = 0, ECurrentState = 1 };
	EState gState;
	
#ifdef _AUTODIFF	
	
	typedef TFad<2, REAL> BFadREAL;	
    
#endif	
	
public:
	
	bool fnewWS;
	
    TPZMultiphase();
    
    TPZMultiphase(int matid, int dim);
    
	virtual ~TPZMultiphase();
    
    /** @brief copy constructor */
    TPZMultiphase(const TPZMultiphase &copy);
    
    TPZMultiphase &operator=(const TPZMultiphase &copy);
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMultiphase"; }
	
	virtual int Dimension() const;
	
	virtual int NStateVariables();	
    
    virtual int MatId();
	
	
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);	
	
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);	
	
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);	
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef);	
	
	virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec,TPZVec<TPZMaterialData> &dataleftvec,TPZVec<TPZMaterialData> &datarightvec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);	
	
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);	
	
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
	
	
	virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
	{
        TPZDiscontinuousGalerkin::Solution(data,dataleftvec,datarightvec,var,Solout,Left,Right);
    }
	
	
	// Here is needed to spefify which models are used in this bi-phasic model see Advanced-Petroleum-Reservoir-Simulation M. Rafiqul Islam
	
	/** @brief K map */	
	TPZStack< TPZFMatrix<REAL> > fKabsoluteMap;	
	
	/** @brief Use or not K map */		
	bool fYorN;
	
	/** @brief Simulation time step */
	REAL fDeltaT;
	
	/** @brief Parameter representing temporal scheme */
	REAL fTheta;	
	
	/** @brief Defines simulation time step. */
	void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}	
	
	/** @brief Defines simulation time step. */
	void SetTheta(REAL timetheta){ this->fTheta = timetheta;}		
	
	void SetLastState(){ gState = ELastState;}
	
	void SetCurrentState(){ gState = ECurrentState;}
	
	/** @brief Defines simulation use a K map */
	void SetYorN(bool dummybool){ this->fYorN = dummybool;}		
	
	
	/**
	 * @name Getting Data
	 * @{
	 */
	
	void SetKMap(TPZStack< TPZFMatrix<REAL> > KabsoluteMap)
	{
		fKabsoluteMap = KabsoluteMap;
	}	
	
	/** @brief Capilar pressure. \f$ pc = pc( Sw ) \f$ */
	void CapillaryPressure(REAL So, REAL &pc, REAL &DpcDSo);
	
	/**
	 * @brief Oil relative permeability.
	 * \f$ Kro = Kro( Sw ) \f$
	 */
	void Kro(REAL Sw, REAL &Kro, REAL &dKroDSw);
	
	/**
	 * @brief Water relative permeability.
	 * \f$ Krw = Krw( Sw ) \f$
	 */
	void Krw(REAL Sw, REAL &Krw, REAL &dKrwSo);
	
	/** 
	 * @brief \f$ Rock porosity. \f$ Phi = Phi( p ) \f$
	 * @param po Refrence pressure
	 */	
	void Porosity(REAL po, REAL &poros, REAL &dPorosDp);
	
	/** 
	 * @brief \f$ Oil density RhoOil = RhoOil( po ) \f$
	 * @param po Refrence pressure
	 */
	void RhoOil(REAL po, REAL &RhoOil, REAL &dRhoOilDpo);
	
	/** 
	 * @brief \f$ Water density RhoWater = RhoWater( pw ) \f$
	 * @param pw Refrence pressure
	 */
	void RhoWater(REAL pw, REAL &RhoWater, REAL &dRhoWaterDpo);
	
	/** 
	 * @brief Oil viscosity. \f$ OilViscosity = ViscOleo( po ) \f$
	 * @param po Refrence pressure
	 */
	void OilViscosity(REAL po, REAL &OilViscosity, REAL &dOilViscosityDpo);
	
	/** 
	 * @brief Water viscosity. \f$ WaterViscosity = WaterViscosity( pw ) \f$
	 * @param po Refrence pressure
	 */	
	void WaterViscosity(REAL po, REAL &WaterViscosity, REAL &dWaterViscosityDpo);
	
	/**
	 * @brief Oil mobility.
	 * \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$
	 */
	void OilLabmda(REAL &OilLabmda, REAL Po, REAL Sw, REAL &dOilLabmdaDPo, REAL &dOilLabmdaDSw);
	
	/**
	 * @brief Water mobility.
	 * \f$ \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$
	 */
	void WaterLabmda(REAL &WaterLabmda, REAL Pw, REAL Sw, REAL &dWaterLabmdaDPw, REAL &dWaterLabmdaDSw);
	
	/**
	 * @brief Bulk mobility.
	 * \f$ \lambda = \lambda( pw , Sw ) \f$
	 */
	void Labmda(REAL &Labmda, REAL Pw, REAL Sw, REAL &dLabmdaDPw, REAL &dLabmdaDSw);
	
	/**
	 * @brief Fractional oil flux.
	 * \f$ f_{Oil} = f_{Oil}( po , Sw ) \f$
	 */
	void fOil(REAL &fOil, REAL Pw, REAL Sw, REAL &dfOilDPw, REAL &dfOilDSw);
	
	/**
	 * @brief Fractional water flux.
	 * \f$ f_{Water} = f_{Water}( pw , Sw ) \f$
	 */
	void fWater(REAL &fWater, REAL Pw, REAL Sw, REAL &dfWaterDPw, REAL &dfWaterDSw);
	
#ifdef _AUTODIFF	
	
	// Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/** @brief Capilar pressure. \f$ pc = pc( Sw ) \f$ */
	void CapillaryPressure(BFadREAL So, BFadREAL &pc);	
	
	/**
	 * @brief Oil relative permeability.
	 * \f$ Kro = Kro( Sw ) \f$
	 */
	void Kro(BFadREAL Sw, BFadREAL &Kro);
	
	/**
	 * @brief Water relative permeability.
	 * \f$ Krw = Krw( Sw ) \f$
	 */
	void Krw(BFadREAL Sw, BFadREAL &Krw);
	
	/** 
	 * @brief \f$ Rock porosity. \f$ Phi = Phi( p ) \f$
	 * @param po Refrence pressure
	 */	
	void Porosity(BFadREAL po, BFadREAL &poros);	
	
	/** 
	 * @brief \f$ Oil density RhoOil = RhoOil( po ) \f$
	 * @param po Refrence pressure
	 */
	void RhoOil(BFadREAL po, BFadREAL &RhoOil);	
	
	/** 
	 * @brief \f$ Water density RhoWater = RhoWater( pw ) \f$
	 * @param pw Refrence pressure
	 */
	void RhoWater(BFadREAL pw, BFadREAL &RhoWater);		
	
	/** 
	 * @brief Oil viscosity. \f$ OilViscosity = ViscOleo( po ) \f$
	 * @param po Refrence pressure
	 */
	void OilViscosity(BFadREAL po, BFadREAL &OilViscosity);	
	
	/** 
	 * @brief Water viscosity. \f$ WaterViscosity = WaterViscosity( pw ) \f$
	 * @param po Refrence pressure
	 */	
	void WaterViscosity(BFadREAL po, BFadREAL &WaterViscosity);	
	
	/**
	 * @brief Oil mobility.
	 * \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$
	 */
	void OilLabmda(BFadREAL OilLabmda, BFadREAL Po, BFadREAL &Sw);	
	
	/**
	 * @brief Water mobility.
	 * \f$ \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$
	 */
	void WaterLabmda(BFadREAL WaterLabmda, BFadREAL Pw, BFadREAL &Sw);	
	
	
	/**
	 * @brief Bulk mobility.
	 * \f$ \lambda = \lambda( pw , Sw ) \f$
	 */
	void Labmda(BFadREAL Labmda, BFadREAL Pw, BFadREAL &Sw);	
	
	/**
	 * @brief Fractional oil flux.
	 * \f$ f_{Oil} = f_{Oil}( po , Sw ) \f$
	 */
	void fOil(BFadREAL fOil, BFadREAL Pw, BFadREAL &Sw);	
	
	
	/**
	 * @brief Fractional water flux.
	 * \f$ f_{Water} = f_{Water}( pw , Sw ) \f$
	 */
	void fWater(BFadREAL fWater, BFadREAL Pw, BFadREAL &Sw);		
	
#endif	
	
	// Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	
	/** @brief Oil density on standard conditions - kg/m3 */
	REAL RhoOilSC();
	
	/** @brief Water density on standard conditions - kg/m3 */
	REAL RhoWaterSC();
	
	/** @brief Gravity */
	TPZFMatrix<REAL> Gravity();
	
	/** @brief Absolute permeability. */
	void K(TPZFMatrix<REAL> &Kab);
	
	/** @brief Absolute permeability inverse. */
	TPZFMatrix<REAL>  Kinv(TPZFMatrix<REAL> &Kab);	
	
	/** @brief Absolute permeability. */
	void LoadKMap(std::string MaptoRead);	
	
	/** @} */
	
};

#endif
