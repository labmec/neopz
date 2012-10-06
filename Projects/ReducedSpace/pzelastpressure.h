/**
 * @file
 * @brief Contains the declaration of the TPZElastPressure class, this material is used to validate the reduce space.
 */

#ifndef PZ_pzelastpressure_h
#define PZ_pzelastpressure_h

#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

#include <iostream>

/**
 * @ingroup material
 * @brief Material to validate the reduce space. 
 * @author Agnaldo Farias
 * @since 8/27/12
 */
/**
 * \li elasticity equation.
 * \f$  div(T(u)) + fxy = 0 \f$ (Eq. 1)
 *
 * \li pressure equation (1d)
 * \f$ *(-wˆ3/visc)div(grad p) + QL= 0 \f$ (Eq. 2) 
 * \li
 */

class TPZElastPressure: public TPZDiscontinuousGalerkin{
    
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
	
	/** @brief Elasticity modulus */
	REAL fE;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
	/** @brief Problem dimension */
	int fDim;
	
	/**
	 * @brief term that multiplies the Laplacian operator, outflow to the poros medio and right side
	 * @note \f$fw \f$ => abertura da fratura 
	 * @note \f$fvisc \f$ => viscosidade do fluido 
	 * @note \f$fQL \f$ => vazao para o meio poroso 
	 * @note \f$fXf \f$ => vetor de carga 
	 */
    REAL fw, fvisc, fQL;
    REAL fXf;
    
	
	/**
	 * @brief Uses plain stress 
     * @note \f$fPlaneStress = 1\f$ => Plain stress state 
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state 
     */
	int fPlaneStress;
	
	REAL fmatId;
    
    /// timestep [s]
	REAL fTimeStep;
    
    /** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	static EState gState;
    
public:
	TPZElastPressure();
	
	TPZElastPressure(int matid, int dim);
	
	virtual ~TPZElastPressure();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZElastPressure"; }
	
	int Dimension() {return fDim;}
	
	virtual int NStateVariables();
	
    
	/** @brief Parameters of pressure: */
	void SetParameters(REAL &xw, REAL &xvisc, REAL &xQL)
	{   
        fw = xw;
        fvisc = xvisc;
        fQL = xQL;
	}
	
	/** 
	 * @brief Set parameters of elastic material:
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 */
	void SetElasticParameters(REAL E, REAL nu,  REAL fx, REAL fy)
	{
		fE = E;
		fnu = nu;
		ff[0] = fx;
		ff[1] = fy;
	}
	
	/**
	 * @brief Set plane problem  
	 * @param planestress \f$ planestress = 1 \f$ indicates use of plain stress or plain strain
	 * planestress = 1 => Plain stress state 
	 * planestress != 1 => Plain Strain state 
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
	
    /// Set the timestep
	void SetTimeStep(REAL delt)
	{
		fTimeStep = delt;
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
    virtual void ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc);
    virtual void ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
	virtual void ApplyDirichlet_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    virtual void ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc);
    virtual void ApplyMixed_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
	
	/** @name Contribute methods
	 * @{
	 */
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param[in] data
	 * @param[in] dataleft Material data from left element
	 * @param[in] dataright Material data from right element
	 * @param[in] weight 
	 * @param[out] ek is the stiffness matrix
	 * @param[out] ef is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
        DebugStop();
    }
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param[in] data 
	 * @param[in] dataleft Material data from left element
	 * @param[in] weight
	 * @param[out] ek is the stiffness matrix
	 * @param[out] ef is the load vector
	 * @param[in] bc is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                       REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
        DebugStop();
    }
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the stiffness matrix
     * @param[out] ef is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {
		DebugStop();
	}
    
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the stiffness matrix
     * @param[out] ef is the load vector
     * @param[in] bc is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
		DebugStop();
	}
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void ContributePressure(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	/** @} */
};

#endif
