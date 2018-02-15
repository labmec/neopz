//
//  TPZMDPMaterial.h
//  PZ
//
//  Created by Agnaldo Farias on 19/09/14.
//
//

#ifndef __PZ__TPZMDPMaterial__
#define __PZ__TPZMDPMaterial__


/**
 * @file
 * @brief Contains the material to run the double projection method (MDP).
 */

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "TPZMatLaplacian.h"

/**
 * @ingroup original operator
 * @brief we have \f$ -fK Laplac(u) = fXf  \f$
 * @brief so that: \f$a(u,v) = grad(u)*grad(v), f(v) = fxf*v\f$
 */
 
 /**
 * @ingroup target operator (inner product)
 * @brief we have \f$ b(u,v) = grad(u)*grad(v) + u*v\f$
 */


class TPZMDPMaterial : public TPZMatLaplacian {
	
	protected :
    
    /** @brief Direction of the convection operator */
    REAL fC[3];
    
public:
	
    
    TPZMDPMaterial(int nummat, int dim);
    
    TPZMDPMaterial(int matid) : TPZMatLaplacian(matid)
    {
        
    }
    
	TPZMDPMaterial();
    
	TPZMDPMaterial(const TPZMDPMaterial &copy);
    
	virtual ~TPZMDPMaterial();
    
	TPZMDPMaterial &operator=(const TPZMDPMaterial &copy);
    
    
	virtual TPZMaterial * NewMaterial(){
		return new TPZMDPMaterial(*this);
	}
    
    
	virtual void Print(std::ostream & out);
    
	virtual std::string Name() { return "TPZMDPMaterial"; }
    
    /**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @brief Contribute to DPG method
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

    /** @} */
    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    virtual int VariableIndex(const std::string &name);
    
	virtual int NSolutionVariables(int var);
    
    void SetConvection(TPZVec<REAL> conv){
        if(conv.size()!=3) DebugStop();
        fC[0]=conv[0];
        fC[1]=conv[1];
        fC[2]=conv[2];
    }
    
public:
    
    
    
    virtual int ClassId() const{
        return TPZMatLaplacianLagrangeID;
    }
    
	virtual void Write(TPZStream &buf, int withclassid) const;
    
	virtual void Read(TPZStream &buf, void *context);
    
};

#endif /* defined(__PZ__TPZMDPMaterial__) */

