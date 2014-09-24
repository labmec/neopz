//
//  pzelasticSest2D.h
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#ifndef __PZ__pzelasticSest2D__
#define __PZ__pzelasticSest2D__

#include <iostream>

/**
 * @file
 * @brief Contains the TPZElasticityMaterial class which implements a two dimensional elastic material in plane stress or strain.
 */

#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzelasmat.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticityMaterialSest2D : public TPZElasticityMaterial {
	
	public :
    
	/** @brief Default constructor */
	TPZElasticityMaterialSest2D();
	/**
	 * @brief Creates an elastic material with:
	 * @param id material id
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	TPZElasticityMaterialSest2D(int id, REAL E, REAL nu, REAL fx, REAL fy, int plainstress = 1);
    
    TPZElasticityMaterialSest2D(int id);
	
	/** @brief Copies the data of one TPZElasticityMaterial object to another */
	TPZElasticityMaterialSest2D(const TPZElasticityMaterial &copy);
	
	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial * NewMaterial() { return new TPZElasticityMaterial(*this);}
	virtual int ClassId() const;

};



#endif /* defined(__PZ__pzelasticSest2D__) */
