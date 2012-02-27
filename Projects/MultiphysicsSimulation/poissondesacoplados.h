/*
 *  LaplacianosDesacoplados.h
 *  PZ
 *
 *  Created by Agnaldo on 10/18/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef POISSONDESACOPLADOSH
#define POISSONDESACOPLADOSH

#include "pzmaterial.h"
#include "pzdiscgal.h"

#include <iostream>


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
/**
 * \f$ fK1 Laplac(u)  = fXf1 ==> Int{fK1 Grad(u).n v}ds - Int{fK1 Grad(u)Grad(v)}dx = Int{fxf1 v}dx   (Eq. 1)  \f$ 
 *
 * \f$ fK2 Laplac(p)  = fXf2  ==> Int{fK2 Grad(p).n v}ds - Int{fK2 Grad(p)Grad(v)}dx = Int{fxf2 v}dx   (Eq. 2) \f$ 
 */


class TwoUncoupledPoisson : public TPZDiscontinuousGalerkin {

protected:
	/** @brief Forcing function value */
	REAL fXf1;
	REAL fXf2;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Coeficient which multiplies the Laplacian operator. */
	REAL fK1;
	REAL fK2;
	
public:
	TwoUncoupledPoisson();
	
	TwoUncoupledPoisson(int matid, int dim);
	
	virtual ~TwoUncoupledPoisson();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TwoUncoupledPoisson"; }
	
	int Dimension() {return fDim;}
	
	virtual int NStateVariables();
	
	void SetParameters(REAL diff1, REAL diff2){
		fK1= diff1;
		fK2 =diff2;
	}
	
	void GetParameters(REAL &diff1, REAL &diff2){
		diff1 = fK1;
		diff2= fK2;
	}
	
	void SetInternalFlux(REAL flux1, REAL flux2)
	{
		fXf1 = flux1;
		fXf2 = flux2;
	}

	//virtual TPZAutoPointer<TPZMaterial> NewMaterial(){
//		return new TwoUncoupledPoisson(*this);
//	}
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	
	//void ContributeInterface(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
//protected:
//	
//	/**
//     * @brief It return a solution to multiphysics simulation.
//     * @param Sol [in] is the solution 
//     * @param DSol [in] is the Gradient
//     * @param axes  [in] is the points of the coordinate axes
//     * @param var [in] number of solution variables. See  NSolutionVariables() method
//     * @param Solout [out] is the solution vector
//     */	
//	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes, int var,TPZVec<REAL> &Solout);

//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
	
	//virtual int IntegrationRuleOrder(TPZVec<int> elPMaxOrder) const;
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	
	
	
	
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef) {
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
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc){
		DebugStop();
	}
	 	
	
};
#endif
