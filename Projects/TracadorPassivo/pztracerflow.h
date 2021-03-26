//
//  pztracerflow.h
//  PZ
//
//  Created by Agnaldo Farias on 16/10/13.
//
//

#ifndef __PZ__pztracerflow__
#define __PZ__pztracerflow__

#include <iostream>

#include "TPZMaterial.h"
#include "pzreal.h"

/**
 * @ingroup material
 * @author Agnaldo Farias
 * @since 16/10/13
 * @brief Contains the TPZTracerFlow class which implements a Tracer Flow Problem 2D with explicit scheme in time.
 *@brief q: flow, p: pressure and s: saturation
 */
/*
 * \f$ poros*ds/dt + div(q*s) = 0  (Eq. 1, transport problem)  \f$
 * \f$ q + (k/visc)*grad(p) = 0  (Eq. 2, elliptic system)  \f$
 * \f$ div(q) = 0  (Eq. 3, elliptic system) \f$
 *
 */


class TPZTracerFlow : public TPZMaterial {
    
protected:
	/** @brief Forcing function value, like source term, to flux-pressure equation */
	REAL fxfPQ;
    
    /** @brief Forcing function value, like source term, to saturation equation */
	REAL fxfS;
	
    /** @brief Permeability of the porous medium. Coeficient which multiplies the gradient operator*/
	REAL fk;
    
	/** @brief Porosity of the porous medium. Coeficient which multiplies the temporal derivative.*/
	REAL fPoros;
    
    /** @brief fluid viscosity*/
	REAL fVisc;
    
    /** @brief convection term (direction velocity) */
	TPZManVector<REAL> fConvDir;
    
    /** @brief Problem dimension */
	int fDim;
    
    //time step
    REAL fTimeStep;
    
    //time of simulation
    REAL fTimeValue;
    
    int fmatId;
    
    bool fPressureEquationFilter;
    
    //Second order Runge-Kutta
    bool fRungeKuttaTwo;
    
	/** @brief State: n ou n+1 */
	enum EState {ELastState = 0, ECurrentState = 1};
    EState gState;
    
	
public:
    
	TPZTracerFlow();
	
	TPZTracerFlow(int matid, int dim);
    
	virtual ~TPZTracerFlow();
    
    /** @brief copy constructor */
    TPZTracerFlow(const TPZTracerFlow &copy);
    
    TPZTracerFlow &operator=(const TPZTracerFlow &copy);
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZTracerFlow"; }
	
	int Dimension() const {return fDim;}
    
    int MatId()
    {
        return fmatId;
    }
	
	virtual int NStateVariables() const {return 1;}
    
    void SetForcesPressure(REAL fxfPQ);
    
    void SetForcesSaturation(REAL fxfS);
    
    void SetPermeability(REAL perm);
    
    void SetViscosity(REAL visc);
    
    void SetPorosity(REAL poros);
    
    void GetPermeability(REAL &perm);
    

	void SetConvectionDirection(TPZVec<REAL> convdir);
	
	void GetConvectionDirection(TPZVec<REAL> &convdir);
    
    void SetLastState(){
        gState = ELastState;
    }
    
	void SetCurrentState(){
        gState = ECurrentState;
    }
    
    void SetPressureEqFilter(){
        fPressureEquationFilter = true;
    }
    void SetFalsePressureEqFilter(){
         fPressureEquationFilter = false;
    }
    
	void SetTimeStep(REAL delt){
		fTimeStep = delt;
	}
    
	void SetTimeValue(REAL TimeValue){
		fTimeValue = TimeValue;
	}
    
    void GetTimeValue(REAL &TimeValue){
		TimeValue = fTimeValue;
	}
    
    void SetTrueRungeKuttaTwo(){
        fRungeKuttaTwo = true;
    }
    
	/** @name Contribute methods
	 * @{
	 */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since March 2, 2013
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
		DebugStop();
	}

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since March 04, 2013
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
		DebugStop();
	}
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param datavec [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since March 04, 2013
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param datavec [in]
	 * @param dataleft [in]
     * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since March 04, 2013
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    
};

#endif /* defined(__PZ__pztracerflow__) */
