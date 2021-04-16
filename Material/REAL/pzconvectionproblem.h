//
//  pzconvectionproblem.h
//  PZ
//
//  Created by Agnaldo Farias on 4/2/13.
//
//

#ifndef __PZ__pzconvectionproblem__
#define __PZ__pzconvectionproblem__


//#include "pzfmatrix.h"
//#include "pzmaterial.h"
//#include "pzvec.h"
#include "TPZMaterial.h"
#include <iostream>


/**
 * @ingroup material
 * @author Agnaldo de Farias
 * @since 4/2/2013
 * @brief Contains the TPZMatConvectionProblem class which implements a convection problem 2D with time dependence
 */
/*
 * \f$ phi*du/dt +  div(ConvDir*u) = f  (Eq. 1)  \f$
 *
 */


class TPZMatConvectionProblem : public TPZMaterial {
    
protected:
	/** @brief Forcing function value */
	REAL fXf;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Coeficient which multiplies the temporal derivative. */
	REAL fRho;
    
    /** @brief convection term (direction velocity) */
	TPZVec<REAL> fConvDir;
    
    //time step
    REAL fTimeStep;
    
    //time of simulation
    REAL fTimeValue;
    
    int fmatId;
    
    //Second order Runge-Kutta
    bool fRungeKuttaTwo;
    
	/** @brief State: n ou n+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
    EState gState;

	
public:
    
	TPZMatConvectionProblem();
	
	TPZMatConvectionProblem(int matid, int dim);
    
	virtual ~TPZMatConvectionProblem();
    
    /** @brief copy constructor */
    TPZMatConvectionProblem(const TPZMatConvectionProblem &copy);
    
    TPZMatConvectionProblem &operator=(const TPZMatConvectionProblem &copy);
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name()  override { return "TPZMatConvectionProblem"; }
	
	int Dimension() const  override {return fDim;}
    
    int MatId()
    {
        return fmatId;
    }
	
	virtual int NStateVariables() const override;
    
    void SetInternalFlux(REAL flux);
	
	void SetParameters(REAL rho, TPZVec<REAL> &convdir);
	
	void GetParameters(REAL &rho, TPZVec<REAL> &convdir);
    
    void SetLastState(){
        gState = ELastState;
    }
    
	void SetCurrentState(){
        gState = ECurrentState;
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
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since March 2, 2013
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since March 04, 2013
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since March 04, 2013
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
     * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since March 04, 2013
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
protected:
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
    public:
virtual int ClassId() const override;

};

#endif /* defined(__PZ__pzconvectionproblem__) */
