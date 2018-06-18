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
#include "pzdiscgal.h"
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
	
    /** @brief permeability tensor. Coeficient which multiplies the gradient operator*/
	TPZFNMatrix<9,REAL> fTensorK;
    
    /** @brief inverse of the permeability tensor.*/
	TPZFNMatrix<9,REAL> fInvK;
    
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
    
    /** @brief Pointer to forcing function, it is the Permeability and its inverse */
    TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;
    
public:
    TPZMixedPoisson();
    
    TPZMixedPoisson(int matid, int dim);
    
    virtual ~TPZMixedPoisson();
    
    TPZMixedPoisson(const TPZMixedPoisson &cp);
    
    TPZMixedPoisson &operator=(const TPZMixedPoisson &copy);
    
    virtual TPZMaterial * NewMaterial(){
        return new TPZMixedPoisson(*this);
    }
    

    virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMixedPoisson"; }
    
    virtual int NStateVariables();
	
	void SetPermeability(REAL perm) {
		fK = perm;
        fTensorK.Zero();
        fInvK.Zero();
        for (int i=0; i<3; i++) {
            fTensorK(i,i) = perm;
            fInvK(i,i) = 1./perm;
        }
	}
    
    //Set the permeability tensor and inverser tensor
   void SetPermeabilityTensor(TPZFMatrix<REAL> K, TPZFMatrix<REAL> invK){
   
//       if(K.Rows() != fDim || K.Cols() != fDim) DebugStop();
//       if(K.Rows()!=invK.Rows() || K.Cols()!=invK.Cols()) DebugStop();
       if(K.Rows() < fDim || invK.Rows() < fDim)
       {
           DebugStop();
       }
       fTensorK = K;
       fInvK = invK;
   }
	
    void SetViscosity(REAL visc) {
		fvisc = visc;
	}
    
	void GetPermeability(REAL &perm) {
		perm = fK;
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
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */	
     virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
        }
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = true;
            }
        }
    }

    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
	
    
    virtual int NEvalErrors() {return 3;}

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
    public:
virtual int ClassId() const;
};

#endif
