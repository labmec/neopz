//
//  pzelasthybrid.h
//  PZ
//
//  Created by Joao on 18/09/12.
//
//

/**
 * @file
 * @brief Contains the TPZElasticityHybridMaterial class which implements a two dimensional elastic material to hybrid method.
 */


#ifndef __PZ__pzelasthybrid__
#define __PZ__pzelasthybrid__

#include <iostream>

#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzelasmat.h" 

/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material to hybrid method. It is derivate from the class TPZElasticityMaterial.
 */
class TPZElasticityHybridMaterial : public TPZElasticityMaterial {
	
	public :
    
	/** @brief Default constructor */
	TPZElasticityHybridMaterial();
	TPZElasticityHybridMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress = 1);
	
	/** @brief Copies the data of one TPZElasticityHybridMaterial object to another */
	TPZElasticityHybridMaterial(const TPZElasticityHybridMaterial &copy);
	
	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial * NewMaterial() { return new TPZElasticityHybridMaterial(*this);}
	
	/** @brief Default destructor */
	virtual ~TPZElasticityHybridMaterial();
	
	/** @brief Returns the model dimension */
	int Dimension() { return 2;}
	
	/** @brief Returns the number of state variables associated with the material */
	virtual  int NStateVariables();
		
	/** @brief Returns the material name*/
	std::string Name() { return "TPZElasticityHybridMaterial"; }
	
    
	/** @name Contribute methods */
	/** @{ */
	
	/** @brief Calculates the element stiffness matrix */
	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
		
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
//	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef);
	
//	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &left, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	    
	
public:
    	
	/** @brief Set PresStress Tensor */
	void SetPreStress(REAL Sigxx, REAL Sigyy, REAL Sigxy);
    
	virtual int ClassId() const;
	

private:
	/** @brief Elasticity modulus */
	REAL fE;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
	/** @brief Forcing vector */
	REAL ff[3];
	
	/** @brief \f$ G = E/2(1-nu) \f$ */
	REAL fEover21PlusNu;
	
	/** @brief \f$ E/(1-nu) \f$ */
	REAL fEover1MinNu2;
	
	/** @brief Pre Stress Tensor - Sigma XX */
	REAL fPreStressXX;
	
	/** @brief Pre Stress Tensor - Sigma YY */
	REAL fPreStressYY;
	
	/** @brief Pre Stress Tensor - Sigma XY */
	REAL fPreStressXY;
	
	/** @brief Uses plain stress */
	int fPlaneStress;
    
    /** @brief indicates which solution should be used for post processing */
    int fPostProcIndex;

};

#endif /* defined(__PZ__pzelasthybrid__) */
