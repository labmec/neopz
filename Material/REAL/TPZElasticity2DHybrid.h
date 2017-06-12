/**
 * @file
 * @brief Contains the TPZElasticityMaterial class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef TPZELASTICITYMATERIAL2DHYBRID
#define TPZELASTICITYMATERIAL2DHYBRID

#include <iostream>

#include "pzmaterial.h"
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
    

	
	/** @brief Applies the element boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    virtual void ContributeVecShapeBC(TPZMaterialData &data,REAL weight,
                              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	
	/** @} */
	
	
public:

    
	virtual int ClassId() const;
	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid);
	
};

#endif
