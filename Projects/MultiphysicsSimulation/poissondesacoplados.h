/**
 * @file
 */

#ifndef POISSONDESACOPLADOSH
#define POISSONDESACOPLADOSH

#include "TPZMaterial.h"


#include <iostream>


/**
 * @ingroup material
 * @author Agnaldo de Farias
 * @since 10/18/2011
 * @brief Contains the TPZMatPoissonDesacoplado class which implements an uncoupled system of the two Poisson equation 
*/
/*
 * \f$ fK1 Laplac(u)  = fXf1 ==> Int{fK1 Grad(u).n v}ds - Int{fK1 Grad(u)Grad(v)}dx = Int{fxf1 v}dx   (Eq. 1)  \f$ 
 *
 * \f$ fK2 Laplac(p)  = fXf2  ==> Int{fK2 Grad(p).n v}ds - Int{fK2 Grad(p)Grad(v)}dx = Int{fxf2 v}dx   (Eq. 2) \f$ 
 */


class TPZMatPoissonDesacoplado : public TPZMaterial {

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
	TPZMatPoissonDesacoplado();
	
	TPZMatPoissonDesacoplado(int matid, int dim);
	
	virtual ~TPZMatPoissonDesacoplado();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatPoissonDesacoplado"; }
	
	int Dimension() const {return fDim;}
	
	virtual int NStateVariables() const {return 1;}
	
	
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

	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	//void ContributeInterface(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */	
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
	
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
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

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
	 	
	
};
#endif
