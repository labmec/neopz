/*
 *  pzporoelastic2d.h
 *  PZ
 *
 *  Created by Agnaldo on 11/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef POISSONDESACOPLADOSH
#define POISSONDESACOPLADOSH

#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"

#include <iostream>


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
/**
 **@ingroup equacao da elasticidade
 * \f$  div(T(u))  + fXf2 = 0  ==> Int{Grad(v).T(u)}dx - Int{v.gN}ds  = Int{ff.v}dx   (Eq. 1) \f$
 *
 *\f$ T(u) =  tr(E(u) − alpha*p*I)lambda*I + 2*nu*(E(u) − alpha*p*I)\f$
 *
 *@ingroup equacao da pressao
 * \f$ -1/visc*div(gradu k)  = 0 ==> k/visc*Int{Grad(u)Grad(v)}dx - Int{k/visv*Grad(u).n v}ds  = 0   (Eq. 2)  \f$ 
 *
 */
 


class TPZPoroElastic2d : public TPZDiscontinuousGalerkin {
	
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
	
	/** @brief Elasticity modulus */
	REAL fE;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
	/** @brief constant poroelastic Biot and K is the drained bulk modulus*/
	REAL falpha; //[1/pa]
	
	/** @brief Problem dimension */
	int fDim;
	
	/**@ permeability of the rock and fluid viscosity*/
	REAL fk; 
	REAL fvisc;
	
	/** @brief Uses plain stress 
	*@ fPlaneStress = 1 => Plain stress state 
	*@ fPlaneStress != 1 => Plain Strain state 
	*/
	int fPlaneStress;
	
	REAL fmatId;
	
public:
	TPZPoroElastic2d();
	
	TPZPoroElastic2d(int matid, int dim);
	
	virtual ~TPZPoroElastic2d();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPoroElastic2d"; }
	
	int Dimension() {return fDim;}
	
	virtual int NStateVariables();

	/** @brief Parameters of rock and fluid: */
	void SetParameters(REAL perm, REAL visc)
	{
		fk = perm;
		fvisc = visc;
	}
	
	/** 
	 * @brief Set parameters of elastic material:
	  * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	void SetParameters(REAL E, REAL nu, REAL falpha, REAL fx, REAL fy)
	{
		fE = E;
		fnu = nu;
		ff[0] = fx;
		ff[1] = fy;
	}
	
	/** @brief Set falpha parameter
	 *@param biot : constant poroelastic Biot [dimensionless]
	 *@param bulk : drained bulk modulus [Pa]
	 */
	void SetAlpha(REAL biot, REAL bulk)
	{
		falpha = biot/bulk;
	}
	
	/** @brief Set plane problem  
	 *@ planestress = 1 => Plain stress state 
	 *@ planestress != 1 => Plain Strain state 
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
	
	//
//	virtual TPZAutoPointer<TPZMaterial> NewMaterial(){
//		return new TPZPoroElastic2d(*this);
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
		
//	virtual int IntegrationRuleOrder(TPZVec<int> elPMaxOrder) const;
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	

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
