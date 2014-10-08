//
//  pzelastoplasticSest2D.h
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#ifndef __PZ__pzelastoplasticSest2D__
#define __PZ__pzelastoplasticSest2D__

#include <iostream>



/**
 * @file
 */


#include "pzmaterial.h"
#include "pzmatwithmem.h"
#include "pzelastoplasticmem.h"
#include "pzporoelastoplasticmem.h"
#include "pzelastoplastic2D.h"
#include "pzmaterial.h"

/**
 * Implements an elastoplastic material and uses the memory feature to store the damage variables
 * This material works only together with the Plasticity Library.
 */

template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZMatElastoPlasticSest2D : public TPZMatElastoPlastic2D<T,TMEM> //, TPZMatWithMem<TMEM>
{
public:
	
	//enum SOLUTIONVARS{ENone = -1};
	/**
	 * Default constructor
	 */
	TPZMatElastoPlasticSest2D();
	
	/** Creates a material object and inserts it in the vector of
	 *  material pointers of the mesh. Upon return vectorindex
	 *  contains the index of the material object within the
	 *  vector
	 */
	TPZMatElastoPlasticSest2D(int id ,  int PlaneStrainOrPlaneStress, STATE sigmaZ);
	
	/** Creates a material object based on the referred object and
	 *  inserts it in the vector of material pointers of the mesh.
	 *  Upon return vectorindex contains the index of the material
	 *  object within the vector
	 */
	TPZMatElastoPlasticSest2D(const TPZMatElastoPlastic2D<T,TMEM> &mat);	
	
    /** @brief Creates a new material from the current object  */
    virtual TPZMaterial * NewMaterial() {
        return new TPZMatElastoPlasticSest2D<T,TMEM>(*this);
    }
    

	/**
	 * Unique identifier for serialization purposes
	 */
    virtual int ClassId() const;
};

#endif /* defined(__PZ__pzelastoplasticSest2D__) */