/**
 * @file
 * @brief Contains the TPZElasticityMaterial class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef TPZELASTICITYMATERIAL2DHYBRID
#define TPZELASTICITYMATERIAL2DHYBRID

#include <iostream>

#include "TPZMaterial.h"
#include "pzelasmat.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticity2DHybrid : public TPZElasticityMaterial {
	
	public :

	/** @brief Default constructor */
	TPZElasticity2DHybrid();
	/** 
	 * @brief Creates an elastic material with:
	 * @param id material id
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	TPZElasticity2DHybrid(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);
    
    TPZElasticity2DHybrid(int id);
	
	/** @brief Copies the data of one TPZElasticityMaterial object to another */
	TPZElasticity2DHybrid(const TPZElasticity2DHybrid &copy);
	
	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial * NewMaterial() { return new TPZElasticity2DHybrid(*this);}
	
	/** @brief Default destructor */
	virtual ~TPZElasticity2DHybrid();
	
    /** @name Contribute methods */
    /** @{ */
    

    virtual void Contribute(TPZVec<TPZMaterialData> &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	/** @brief Applies the element boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
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
  //  virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    

	/** @} */
	
    /** @name Post processing methods
     * @{
     */
    

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
  //  virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** @} */
    
//    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
//    {
//        TPZDiscontinuousGalerkin::Errors(data[0],u_exact,du_exact,errors);
//    }

public:

    
	public:
virtual int ClassId() const;

	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid) const;
	
};

#endif
