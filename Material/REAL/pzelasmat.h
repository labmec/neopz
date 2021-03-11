/**
 * @file
 * @brief Contains the TPZElasticityMaterial class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef ELASMATHPP
#define ELASMATHPP

#include <iostream>

#include "TPZMaterial.h"



/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticityMaterial : public TPZMaterial {
	
	public :

	/** @brief Default constructor */
	TPZElasticityMaterial();
	/** 
	 * @brief Creates an elastic material with:
	 * @param id material id
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	TPZElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);
    
    TPZElasticityMaterial(int id);
	
	/** @brief Copies the data of one TPZElasticityMaterial object to another */
	TPZElasticityMaterial(const TPZElasticityMaterial &copy);
	
	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial * NewMaterial()  override { return new TPZElasticityMaterial(*this);}
	
	/** @brief Default destructor */
	virtual ~TPZElasticityMaterial();
	
    /** @brief Set elasticity parameters */
    void SetElasticity(REAL E, REAL nu)
    {
        fE_def	= E;  // Young modulus
        fnu_def	= nu;   // poisson coefficient

    }
    
    /// Set a variable elasticity and poisson coefficient
    void SetElasticityFunction(TPZAutoPointer<TPZFunction<STATE> > func)
    {
        fElasticity = func;
    }
    
    /// Set the material configuration to plane strain
    void SetPlaneStrain()
    {
        fPlaneStress = 0;
    }
    
    /// Set the material configuration to plane stress
    void SetPlaneStress()
    {
        fPlaneStress = 1;
    }
    
    /** @brief Set forcing function */
    void SetBodyForce(REAL fx, REAL fy)
    {
        ff[0] = fx;
        ff[1] = fy;
        ff[2] = 0.;
    }
	/** @brief Returns the model dimension */
	int Dimension() const  override { return 2;}
	
	/** @brief Returns the number of state variables associated with the material */
	virtual  int NStateVariables() const override;
	
	/** @brief Print the material data*/
	virtual void Print(std::ostream & out = std::cout) override;
	
	/** @brief Returns the material name*/
	std::string Name()  override { return "TPZElasticityMaterial"; }
	
	/** @brief Returns the number of components which form the flux function */
	virtual short NumberOfFluxes(){return 3;}
		
	/** @name Contribute methods */
	/** @{ */
	
	/** @brief Calculates the element stiffness matrix */
	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
//    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ef)
//    {
//        DebugStop();
//    }
    
    void ContributeVecShape(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	/** @brief Calculates the element stiffness matrix */
	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;

    /** @brief Applies the element boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    void ContributeVecShapeBC(TPZMaterialData &data,REAL weight,
                              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	/** @brief Applies the element boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
    
    //virtual void FillDataRequirements(TPZMaterialData &data);
     virtual void FillDataRequirements(TPZMaterialData &data) override;
    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data) override;
        
     
    
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
		PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override {
		PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &left, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	}
	
	/** @} */
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name) override;
	
	/** 
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 */
	virtual int NSolutionVariables(int var) override;
    
    STATE GetLambda(REAL E, REAL nu) const
    {
        STATE lambda = (nu*E)/((1.+nu)*(1.-2.*nu));
        return lambda;
    }
    
    STATE GetMU(REAL E, REAL nu) const
    {
        STATE mu = E/(2.*(1.+nu));
        return mu;
    }
	
public:

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::SolutionDisc(data,dataleft,dataright,var,Solout);
	}
	
    /** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
    virtual int NEvalErrors() override {return 6;}
    

    
	/** 
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
	virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;//Cedric
	
	/** @brief Returns the elasticity modulus E */
	REAL E() {return fE_def;}
	
	/** @brief Returns the poison coefficient modulus E */
	REAL Nu() {return fnu_def;}
	
	/** @brief Set PresStress Tensor */
	void SetPreStress(REAL Sigxx, REAL Sigyy, REAL Sigxy, REAL Sigzz);
    
	public:
virtual int ClassId() const override;

	
	virtual void Read(TPZStream &buf, void *context) override;
	
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	
	
protected:
	/** @brief Elasticity modulus */
	REAL fE_def;
	
	/** @brief Poison coeficient */
	REAL fnu_def;
	
    /** Elasticity function */
    TPZAutoPointer<TPZFunction<STATE> > fElasticity;
    
	/** @brief Forcing vector */
	TPZManVector<STATE,3> ff;
	
	/** @brief \f$ G = E/2(1-nu) \f$ */
	REAL fEover21PlusNu_def;
	
	/** @brief \f$ E/(1-nu) \f$ */
	REAL fEover1MinNu2_def;
	
	/** @brief Pre Stress Tensor - Sigma XX */
	REAL fPreStressXX;
	
	/** @brief Pre Stress Tensor - Sigma YY */
	REAL fPreStressYY;
	
	/** @brief Pre Stress Tensor - Sigma XY */
	REAL fPreStressXY;
    
    /** @brief Pre Stress Tensor - Sigma ZZ */
    REAL fPreStressZZ;
	
	/** @brief Uses plain stress */
	int fPlaneStress;
    
};

#endif
