//
//  mixedpoisson.h
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_mixedpoisson_h
#define PZ_mixedpoisson_h

#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzpoisson3d.h"
#include "pzmaterial.h"
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
	
	/** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
	REAL fk;
    
    /** @brief permeability tensor. Coeficient which multiplies the gradient operator*/
	TPZFMatrix<REAL> fTensorK;
    
    /** @brief inverse of the permeability tensor.*/
	TPZFMatrix<REAL> fInvK;
    
    /** @brief fluid viscosity*/
	REAL fvisc;
    
    /** @brief Problem dimension */
	int fDim;
    
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
    
    
    virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMixedPoisson"; }
    
    virtual int NStateVariables();
	
	void SetPermeability(REAL perm) {
		fk = perm;
	}
    
    //Set the permeability tensor and inverser tensor
   void SetPermeabilityTensor(TPZFMatrix<REAL> K, TPZFMatrix<REAL> invK){
   
       if(K.Rows() != fDim || K.Cols() != fDim) DebugStop();
       if(K.Rows()!=invK.Rows() || K.Cols()!=invK.Cols()) DebugStop();
       
       fTensorK = K;
       fInvK = invK;
   }
	
    void SetViscosity(REAL visc) {
		fvisc = visc;
	}
    
	void GetPermeability(REAL &perm) {
		perm = fk;
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
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
	
    	
};

#endif
