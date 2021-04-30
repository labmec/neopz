//
//  mixedpoisson.h
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef TPZMIXEDPOISSONPARABOLICH
#define TPZMIXEDPOISSONPARABOLICH

#include "mixedpoisson.h"

/**
 * @ingroup material
 * @author Philippe Devloo
 * @since Sept 4, 2017
 * @brief Material to solve a mixed time dependent poisson problem 2d by multiphysics simulation
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
/**
 * \f$ Q = -(k/visc)*grad(p)  ==> Int{Q.q}dx - (k/visc)*Int{p*div(q)}dx + (k/visc)*Int{pD(q.n)}ds = 0  (Eq. 1)  \f$ 
 *
 * \alpha p^{n+1)/\Delta t div(Q) = f +\alpha p^n/\Delta t ==>
 * \f$ Int{\alpha v*p^{n+1} +div(Q)*v}dx = Int{v*p^n*\alpha /\Delta t f^{n+1}*v}dx (Eq. 2) \f$
 *
 * \f$ p = pD in Dirichlet boundary and Q.n = qN in Neumann boundary\f$
 */


class TPZMixedPoissonParabolic : public TPZMixedPoisson {
    
protected:
    
    // thermal capacity or compressibility
    STATE fRho;
    
    // timestep used
    STATE fDeltaT;
    
public:
    TPZMixedPoissonParabolic();
    
    TPZMixedPoissonParabolic(int matid, int dim);
    
    virtual ~TPZMixedPoissonParabolic();
    
    TPZMixedPoissonParabolic(const TPZMixedPoissonParabolic &cp);
    
    TPZMixedPoissonParabolic &operator=(const TPZMixedPoissonParabolic &copy);
    
    void SetParabolicConstants(STATE rho, STATE deltat)
    {
        fRho = rho;
        fDeltaT = deltat;
    }
    
    void SetDeltaT(STATE deltat)
    {
        fDeltaT = deltat;
    }
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
        if (type == 3) {
            data.fNeedsNormal = true;
        }
    }

    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
        }
        datavec[1].fNeedsSol = true;
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = true;
            }
        }
    }

    
    virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override { return "TPZMixedPoissonParabolic"; }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

    
};

#endif
