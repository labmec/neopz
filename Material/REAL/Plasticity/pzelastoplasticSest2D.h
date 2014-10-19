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
	TPZMatElastoPlasticSest2D(const TPZMatElastoPlasticSest2D<T,TMEM> &mat);
	
    /** @brief Creates a new material from the current object  */
    virtual TPZMaterial * NewMaterial() {
        return new TPZMatElastoPlasticSest2D<T,TMEM>(*this);
    }
    

	/**
	 * Unique identifier for serialization purposes
	 */
    virtual int ClassId() const;
    
    /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
     * @param DeltaStrain [out]
     * @param data [in]
     */
    void ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain);
    
    /** Calls the plasticity template aggregate applyStrainComputeDep method
     *  @param data [in]
     *  @param DeltaStrain [in]
     *  @param Stress [out]
     *  @param Dep [out]
     */
    void ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
                                    TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);
    

        
    /**
     * Set the deformation in the Z direction
     */
    void SetZDeformation(STATE epsz)
    {
        fZDeformation = epsz;
    }
    
    /**
     * @brief return the increment in z deformation
     */
    STATE ZDeformation()
    {
        return fZDeformation;
    }

private:
    /**
     * @brief ZDeformation of the mesh
     */
    STATE fZDeformation;
    


};

#endif /* defined(__PZ__pzelastoplasticSest2D__) */