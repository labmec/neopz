/**
 * \file
 * @brief Contains the TPZL2Projection class which implements an L2 projection to constant solution values.
 */
//$Id: pzl2projection.h,v 1.14 2011-05-26 03:28:57 phil Exp $

#ifndef PZL2PROJECTION_H
#define PZL2PROJECTION_H

#include "pzmaterial.h"
#include "pzdiscgal.h"

/**
 * @ingroup material
 * @brief Implements an L2 projection to constant solution values.
 * @since April 23, 2007
 */
class TPZL2Projection : public TPZDiscontinuousGalerkin{
	
protected:
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Number of state variables */
	int fNStateVars;
	
	/** @brief Constant solution vector */
	TPZVec<REAL> fSol;
	
	/** @brief Argument defining this material is a referred material */
	bool fIsReferred;
	
	/** @brief Order for setting the integration rule */
	int fIntegrationOrder;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale;
	
public:
	
	/**
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param sol constant solution vector
	 * @param IntegrationOrder numeric integration order
	 */
	TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol,
					int IntegrationOrder = -1);
	
	/** @brief Default destructor */
	~TPZL2Projection();
	
	/** @brief Copy constructor */
	TPZL2Projection(const TPZL2Projection &cp);
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 };
	
    /** 
     * @brief Get the order of the integration rule necessary to integrate an
     * element with polinomial order p
     */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const;

    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {
        fScale = scale;
    }
    
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef)
	{
		TPZDiscontinuousGalerkin::Contribute(data,weight,ef);
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ef){
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ef,TPZBndCond &bc){
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief Returns problem dimension */
	virtual int Dimension(){ return this->fDim; }
	
	/** @brief Returns number of state variables */
	virtual int NStateVariables(){ return this->fNStateVars; }
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ef, TPZBndCond &bc)
	{
		TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
	}
	
	/** @brief Define if material is referred or not */
	void SetIsReferred(bool val);
	
	/** @brief To create another material of the same type */
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	/** @brief It returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
protected:
	/** @brief It returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
						  TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation*/
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}

	
    /** @brief Returns the solution associated with the var index based on the finite element approximation*/
	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<REAL> &Solout)
	{
		TPZDiscontinuousGalerkin::SolutionDisc(data,dataleft,dataright,var,Solout);
	}

};

#endif
