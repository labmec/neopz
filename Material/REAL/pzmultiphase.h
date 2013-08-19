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

/**
 * @ingroup material
 * @author Omar Duran
 * @since 19/08/2013
 * @brief Material to solve a 2d multiphase transport problems by multiphysics simulation
 * @brief Here is used H1, Hdiv ... spaces
 */


class TPZMultiphase : public TPZDiscontinuousGalerkin {
    
protected:
	/** @brief Definition of constants */
	REAL ff;
    
public:
    TPZMultiphase();
    
    TPZMultiphase(int matid, int dim);
    
	virtual ~TPZMultiphase();
    
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMultiphase"; }
    
	virtual int NStateVariables();
	
    
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
	
    /**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */	
	
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
	    	
};

#endif
