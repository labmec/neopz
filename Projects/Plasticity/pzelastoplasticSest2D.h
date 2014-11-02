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
	  EStrainXX = 1,
	  EStrainYY = 2,
	  EStrainZZ = 3,
	  EStrainXY = 4,
	  EStrainXZ = 5,
	  EStrainYZ = 6,
	  // Elastic Strain
	  EElStrainVol = 7,
	  EElStrainXX = 8,
	  EElStrainYY = 9,
	  EElStrainZZ = 10,
	  EElStrainXY = 11,
	  EElStrainXZ = 12,
	  EElStrainYZ = 13,
	  // Plastic Strain
	  EPlStrainVol = 14,
	  EPlStrainXX = 15,
	  EPlStrainYY = 16,
	  EPlStrainZZ = 17,
	  EPlStrainXY = 18,
	  EPlStrainXZ = 19,
	  EPlStrainYZ = 20,
	  EPlStrainSqJ2 = 21,
	  EPlStrainSqJ2El = 22,
	  EPlAlpha = 23,
	  // Displacement
	  EDisplacementX = 24,
	  EDisplacementY = 25,
	  EDisplacementZ = 26,
	  EDisplacementTotal = 27,
	  // Total Stress
	  ETotStressI1 = 28,
	  ETotStressJ2 = 29,
	  ETotStressXX = 30,
	  ETotStressYY = 31,
	  ETotStressZZ = 32,
	  ETotStressXY = 33,
	  ETotStressXZ = 34,
	  ETotStressYZ = 35,
	  ETotStress1 = 36,
	  ETotStress2 = 37,
	  ETotStress3 = 38,
	  // Effective stress
	  EEffStressI1 = 39,
	  EEffStressJ2 = 40,
	  EEffStressXX = 41,
	  EEffStressYY = 42,
	  EEffStressZZ = 43,
	  EEffStressXY = 44,
	  EEffStressXZ = 45,
	  EEffStressYZ = 46,
	  EEffStress1 = 47,
	  EEffStress2 = 48,
	  EEffStress3 = 49,
	  // Yield Surface
	  EYieldSurface1 = 50,
	  EYieldSurface2 = 51,
	  EYieldSurface3 = 52,
	  // Simulation
	  EPOrder = 53,
	  ENSteps = 54,
	  // Pore pressure
	  EPorePressure = 55,
	  // Material
	  EMatPorosity = 56,
	  EMatE = 57,
	  EMatPoisson = 58
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
	TPZMatElastoPlasticSest2D(int id ,  int PlaneStrainOrPlaneStress);
	
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
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    virtual void Read(TPZStream &buf, void *context);
    

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

    /**
     * Set the biot coefficient
     */
    void SetBiot(STATE biot)
    {
      fbiot = biot;
    }
    
    /**
     * @brief return the biot coefficient
     */
    STATE Biot()
    {
      return fbiot;
    }
  
    /**returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name);
    
      /** returns the number of variables associated with the variable
	  indexed by var.  var is obtained by calling VariableIndex*/
      virtual int NSolutionVariables(int var);

    
    
    /**returns the solution associated with the var index based on
     * the finite element approximation*/
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);

private:
    /**
     * @brief ZDeformation of the mesh
     */
    STATE fZDeformation;
  
    /// constante de Biot
    STATE fbiot;


};

#endif /* defined(__PZ__pzelastoplasticSest2D__) */