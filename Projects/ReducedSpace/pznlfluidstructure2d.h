//
//  pznlfluidstructure2d.h
//  PZ
//
//  Created by Agnaldo Farias on 9/17/12.
//
//

#ifndef __PZ__pznlfluidstructure2d__
#define __PZ__pznlfluidstructure2d__

#include <iostream>

#include "TPZMaterial.h"

#include "pzvec.h"
#include "pzfmatrix.h"

#include "pznlfluidstructureData.h"

#include <iostream>

/**
 * @ingroup material
 * @brief Material  to validate the reduce space.
 */
/**
 **@ingroup elasticity equation.
 * \f$  div(T(u)) + fxy = 0 \f$ (Eq. 1)
 *
 *@ingroup pressure equation (1d)
 * \f$ *(-w^3/visc)div(grad p) + QL= 0 (Eq. 2)  \f$
 *
 */


class TPZNLFluidStructure2d : public TPZMaterial
{
    
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
    
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief term that multiplies the Laplacian operator, outflow to the poros medio and right side
     * @note \f$fXf f$ => vetor de carga
     */
    REAL fXf;
	
	/** @brief Uses plain stress
     * @note \f$fPlaneStress = 1\f$ => Plain stress state
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state
     */
	int fPlaneStress;
	
	REAL fmatId;
    REAL fE;//young
    REAL fPoiss;//poisson
    REAL fVisc;//viscosity
    
    /** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	static EState gState;
    
public:
    
	TPZNLFluidStructure2d();
	
	TPZNLFluidStructure2d(int matid, int dim, REAL young, REAL poiss, REAL visc);
	
	virtual ~TPZNLFluidStructure2d();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZNLFluidStructure2d"; }
	
	int Dimension() const {return fDim;}
	
	virtual int NStateVariables() const {return 1;}
	
	/** @brief Set plane problem
	 * planestress = 1 => Plain stress state
	 * planestress != 1 => Plain Strain state
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
    
    int MatId()
    {
        return fmatId;
    }
    
    
    void SetLastState(){ gState = ELastState; }
	void SetCurrentState(){ gState = ECurrentState; }
    
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
    
    virtual void ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    virtual void ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc);
    virtual void ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    virtual void ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
    {
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
                                       REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
    {
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
    {
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
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc)
    {
		DebugStop();
	}
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    void ContributePressure(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
};

#endif /* defined(__PZ__pznlfluidstructure2d__) */
