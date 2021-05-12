/**
 * @file
 * @brief Contains the TPZL2Projection class which implements an L2 projection to constant solution values.
 */

#ifndef PZL2PROJECTION_H
#define PZL2PROJECTION_H

#include "TPZMaterial.h"


/**
 * @ingroup material
 * @brief Implements an L2 projection to constant solution values.
 * @since April 23, 2007
 */
class TPZL2Projection : public TPZMaterial{
	
protected:
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Number of state variables */
	int fNStateVars;
	
	/** @brief Constant solution vector */
	TPZVec<STATE> fSol;
	
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
	TPZL2Projection(int id, int dim, int nstate, TPZVec<STATE> &sol,
					int IntegrationOrder = -1);
	
	/** @brief Default destructor */
	~TPZL2Projection();
	
	/** @brief Copy constructor */
	TPZL2Projection(const TPZL2Projection &cp);
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , EDerivative = 2};
	
    /** 
     * @brief Get the order of the integration rule necessary to integrate an
     * element with polinomial order p
     */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const override ;
	
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {
        fScale = scale;
    }
    
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override ;
    virtual void ContributeVecShape(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override {
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		//NOTHING TO BE DONE HERE
	}
	
	/** @brief Returns problem dimension */
	virtual int Dimension() const  override { return this->fDim; }
    
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }
	
	/** @brief Returns number of state variables */
	virtual int NStateVariables() const override { return this->fNStateVars; }
    
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override ;
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {
		TPZMaterial::ContributeBC(data,weight,ef,bc) ;
	}
	
	/** @brief Define if material is referred or not */
	void SetIsReferred(bool val);
	
	/** @brief To create another material of the same type */
	virtual TPZMaterial * NewMaterial() override ;
	
	/** @brief It returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name) override ;
	
	virtual int NSolutionVariables(int var) override ;
	
protected:
	/** @brief It returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol,
						  TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout) override ;
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation*/
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override 
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
    /** @brief Returns the solution associated with the var index based on the finite element approximation*/
	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::SolutionDisc(data,dataleft,dataright,var,Solout);
	}
protected:
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
           TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
           TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &    values) override;
    public:
virtual int ClassId() const override ;

};

#endif
