
/**
 * @file
 * @brief Contains the TPZL2Projection class which implements an L2 projection to use in reconstruction gradient.
 */

//
//  pznewl2projection.h
//  PZ
//
//  Created by Agnaldo Farias on 1/9/13.
//
//



#ifndef PZ_pznewl2projection_h
#define PZ_pznewl2projection_h

#include "pzmaterial.h"
#include "pzdiscgal.h"

class TPZL2ProjectionForGradient : public TPZDiscontinuousGalerkin{

protected:
    
    /** @brief Problem dimension */
	int fDim;
	
	/** @brief Number of state variables */
	int fNStateVars;
	
    REAL fmatId;
    
	/**
     *@brief data of gradient vector:
     * first column: index of the element
     * second column: midpoint of the element (x0, y0, z0)
     * other columns: components of the gradient reconstructed (dudx, dudy, dudz)
     */
	TPZFMatrix<REAL> fgradients;
    
    public:
    
    /**
	 * @brief Class constructor
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param gradients data of gradient vector
	 */
	TPZL2ProjectionForGradient(int matid, int dim,int nstate);
	
	/** @brief Default destructor */
	~TPZL2ProjectionForGradient();
	    
	
	void SetGradients(TPZFMatrix<REAL> &grad) {
		fgradients = grad;
	}
	
    /** @brief Returns problem dimension */
	virtual int Dimension(){ return this->fDim; }
	
	/** @brief Returns number of state variables */
	virtual int NStateVariables(){ return this->fNStateVars; }
    
    int MatId()
    {
        return fmatId;
    }
    
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
    
	//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
        DebugStop();
    }
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                       REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
        DebugStop();
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
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {
		DebugStop();
	}
    
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
		DebugStop();
	}
};
#endif
