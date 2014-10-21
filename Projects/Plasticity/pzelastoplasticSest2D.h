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
  
        enum SOLUTIONVARS{ENone = -1,
	  // Strain
	  EStrainVol = 0,
	  EStrainXX,
	  EStrainYY,
	  EStrainZZ,
	  EStrainXY,
	  EStrainXZ,
	  EStrainYZ,
	  // Elastic Strain
	  EElStrainVol,
	  EElStrainXX,
	  EElStrainYY,
	  EElStrainZZ,
	  EElStrainXY,
	  EElStrainXZ,
	  EElStrainYZ,
	  // Plastic Strain
	  EPlStrainVol,
	  EPlStrainXX,
	  EPlStrainYY,
	  EPlStrainZZ,
	  EPlStrainXY,
	  EPlStrainXZ,
	  EPlStrainYZ,
	  EPlStrainSqJ2,
	  EPlStrainSqJ2El,
	  EPlAlpha,
	  // Displacement
	  EDisplacementX,
	  EDisplacementY,
	  EDisplacementZ,
	  EDisplacementTotal,
	  // Total Stress
	  ETotStressI1,
	  ETotStressJ2,
	  ETotStressXX,
	  ETotStressYY,
	  ETotStressZZ,
	  ETotStressXY,
	  ETotStressXZ,
	  ETotStressYZ,
	  ETotStress1,
	  ETotStress2,
	  ETotStress3,
	  // Effective stress
	  EEffStressI1,
	  EEffStressJ2,
	  EEffStressXX,
	  EEffStressYY,
	  EEffStressZZ,
	  EEffStressXY,
	  EEffStressXZ,
	  EEffStressYZ,
	  EEffStress1,
	  EEffStress2,
	  EEffStress3,
	  // Yield Surface
	  EYieldSurface1,
	  EYieldSurface2,
	  EYieldSurface3,
	  // Simulation
	  EPOrder,
	  ENSteps,	  
	  // Pore pressure
	  EPorePressure,
	  // Material
	  EMatPorosity,
	  EMatE,
	  EMatPoisson
};
	
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
    
    /**returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name);
    
      /** returns the number of variables associated with the variable
	  indexed by var.  var is obtained by calling VariableIndex*/
      virtual int NSolutionVariables(int var);

private:
    /**
     * @brief ZDeformation of the mesh
     */
    STATE fZDeformation;
    


};

#endif /* defined(__PZ__pzelastoplasticSest2D__) */