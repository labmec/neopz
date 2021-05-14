//
//  mixedpoisson.h
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_mixedpoisson_h
#define PZ_mixedpoisson_h

#include "TPZMaterial.h"

#include "pzpoisson3d.h"
#include "TPZMaterial.h"
#include "pzfunction.h"

/**
 * @ingroup material
 * @author Agnaldo Farias
 * @since 5/28/2012
 * @brief Material to solve a mixed poisson problem 2d by multiphysics simulation
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
/**
 * \f$ Q = -(k/visc)*grad(p)  ==> Int{Q.q}dx - (k/visc)*Int{p*div(q)}dx + (k/visc)*Int{pD(q.n)}ds = 0  (Eq. 1)  \f$ 
 *
 * \f$ div(Q) = f  ==> Int{div(Q)*v}dx = Int{f*v}dx (Eq. 2) \f$ 
 *
 * \f$ p = pD in Dirichlet boundary and Q.n = qN in Neumann boundary\f$
 */


class TPZMixedPoisson : public TPZMatPoisson3d {
    
protected:
	/** @brief Forcing function value */
	REAL ff;
    
    /** @brief fluid viscosity*/
	REAL fvisc;
    
    /** @brief Choose Stabilized method */
    bool fIsStabilized;
    
    /** @brief Coeficient of Stabilization*/
    REAL fdelta1;
    REAL fdelta2;
    
    /** @brief Coeficient that multiplies the Stabilization term fdelta1*/
    REAL fh2;
    bool fUseHdois;
    /** @brief post-processing procedure for error estimation as Ainsworth*/
    
    
    /** @brief Pointer to forcing function, it is the Permeability and its inverse */
    TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;
    
public:
    TPZMixedPoisson();
    
    TPZMixedPoisson(int matid, int dim);
    
    virtual ~TPZMixedPoisson();
    
    TPZMixedPoisson(const TPZMixedPoisson &cp);
    
    TPZMixedPoisson &operator=(const TPZMixedPoisson &copy);
    
    virtual TPZMaterial * NewMaterial() override{
        return new TPZMixedPoisson(*this);
    }
    

    virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override{ return "TPZMixedPoisson"; }
    
    virtual int NStateVariables() const override;

    void SetViscosity(REAL visc) {
		fvisc = visc;
	}
    
    void GetMaxPermeability(REAL &perm)
    {
        if(fPermeabilityFunction) DebugStop();//Not implemented for permeability given by a function
        TPZFNMatrix<9,STATE> TensorK(3,3);
        this->GetPermeability(TensorK);
        perm = 0;
        for(int i=0; i<3; i++) perm = perm < TensorK(i,i) ? TensorK(i,i) : perm;
    }
	
	void SetInternalFlux(REAL flux) {
		ff = flux;
	}
    
    void SetStabilizedMethod(){
        fIsStabilized = true;
    }
    
    void SetHdois(){
        fUseHdois = true;
    }
    
    void SetStabilizationCoeficients(REAL delta1, REAL delta2){
        fdelta1 = delta1;
        fdelta2 = delta2;
    }
    
    void SetPermeabilityFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fPermeabilityFunction = fp;
    }

    TPZAutoPointer<TPZFunction<STATE> > PermeabilityFunction()
    {
        return fPermeabilityFunction;
    }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */	
     virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsSol = false;
        }
        datavec[0].fNeedsNormal = true;
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = false;
            }
        }
    }

    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
	
    
    virtual int NEvalErrors() override {return 5;}

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors) override;
    
    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc) override;

    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors,TPZBndCond &bc) override;
    
    
    public:
virtual int ClassId() const  override;

};

#endif
