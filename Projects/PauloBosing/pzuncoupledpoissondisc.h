//
//  pzuncoupledpoissondisc.h
//  PZ
//
//  Created by Agnaldo Farias on 2/19/13.
//
//

#ifndef __PZ__pzuncoupledpoissondisc__
#define __PZ__pzuncoupledpoissondisc__


#include "pzfmatrix.h"

#include <iostream>


/**
 * @ingroup material
 * @author Agnaldo de Farias
 * @since 10/18/2011
 * @brief Contains the TPZMatUncoupledPoissonDisc class which implements an uncoupled system of the two Poisson equation
 */
/*
 * \f$ - fK1 Laplac(u1)  = fXf1 ==> Int{fK1 Grad(u1)Grad(v)}dx  - Int{fK1 Grad(u1).n v}ds = Int{fxf1 v}dx   (Eq. 1)  \f$
 *
 * \f$ - fK2 Laplac(u2)  = fXf2 ==> Int{fK2 Grad(u2)Grad(v)}dx - Int{fK2 Grad(u2).n v}ds  = Int{fxf2 v}dx   (Eq. 2) \f$
 */


class TPZMatUncoupledPoissonDisc : public TPZMaterial {
    
protected:
	/** @brief Forcing function value */
	REAL fXf1;
	REAL fXf2;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Coeficient which multiplies the Laplacian operator. */
	REAL fK1;
	REAL fK2;
    
    /** @brief Symmetry coefficient of elliptic term */
	/**
	 * Symmetrical formulation - Global element method - has coefficient = -1. \n
	 * Non-symmetrical formulation - Baumann's formulation - has coefficient = +1.
	 */
	REAL fSymmetry1;
    REAL fSymmetry2;
    
    /** @brief Constant multiplyer of penalty term, when required is set. */
	REAL fPenaltyConstant1;
    REAL fPenaltyConstant2;
	
public:
    
	TPZMatUncoupledPoissonDisc();
	
	TPZMatUncoupledPoissonDisc(int matid, int dim);
    
	virtual ~TPZMatUncoupledPoissonDisc();
    
    /** @brief copy constructor */
    TPZMatUncoupledPoissonDisc(const TPZMatUncoupledPoissonDisc &copy);
    
    TPZMatUncoupledPoissonDisc &operator=(const TPZMatUncoupledPoissonDisc &copy);
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatUncoupledPoissonDisc"; }
	
	int Dimension() const {return fDim;}
	
	virtual int NStateVariables() const { return 1;}
	
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
    
    /** @brief Set material elliptic term as the global element method, i.e. the symmetrical formulation */
	void SetSymmetryOne(){
		this->fSymmetry1 = -1.0;
	}
	
	/** @brief Set material elliptic term as the Baumann's formulation, i.e. the non-symmetrical formulation */
	void SetNonSymmetricOne() {
		this->fSymmetry1 = +1.0;
	}
    
    /** @brief Set material elliptic term as the global element method, i.e. the symmetrical formulation */
	void SetSymmetryTwo(){
		this->fSymmetry2 = -1.0;
	}
	
	/** @brief Set material elliptic term as the Baumann's formulation, i.e. the non-symmetrical formulation */
	void SetNonSymmetricTwo() {
		this->fSymmetry2 = +1.0;
	}
    
    void SetPenaltyConstant(REAL penalty_const1, REAL penalty_const2) {
		this->fPenaltyConstant1 = penalty_const1;
        this->fPenaltyConstant2 = penalty_const2;
	}
    
    bool IsSymetricOne(){
		if (fSymmetry1 == -1.0) return true;
		if (fSymmetry1 == +1.0) return false;
		PZError << __PRETTY_FUNCTION__ << "\n Comparacao de numeros reais da errado\n";
		return false;
	}
    
    bool IsSymetricTwo(){
		if (fSymmetry2 == -1.0) return true;
		if (fSymmetry2 == +1.0) return false;
		PZError << __PRETTY_FUNCTION__ << "\n Comparacao de numeros reais da errado\n";
		return false;
	}
    
    virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
	
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
		DebugStop();
	}
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
		DebugStop();
	}
    
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief Computes a contribution to the stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since February 21, 2013
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @} */

};

#endif /* defined(__PZ__pzuncoupledpoissondisc__) */
