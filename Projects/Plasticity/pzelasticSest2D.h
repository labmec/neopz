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
#include "TPZMatElasticity2D.h"
#include "TPZTensor.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticityMaterialSest2D : public TPZMatElasticity2D {
    
    
    public :
    
    enum SOLUTIONVARS{
        ENone = -1,
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
        EMatPoisson = 58 };
    
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
    TPZElasticityMaterialSest2D(const TPZElasticityMaterialSest2D &copy);
    
    /// Destructor
    virtual ~TPZElasticityMaterialSest2D();
    
    TPZElasticityMaterialSest2D &operator=(const TPZElasticityMaterialSest2D &copy);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    /** @brief Creates a new material from the current object  */
    virtual TPZMaterial * NewMaterial() { return new TPZElasticityMaterialSest2D(*this);}
    
    virtual int ClassId() const { return TPZElasticityMaterialSest2DID; }
    
    /**
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    virtual void Read(TPZStream &buf, void *context);
    
    
    /** @brief Set M. Biot alpha constant */
    void SetBiotAlpha(REAL Alpha)
    {
        fBiotAlpha = Alpha;
    }
    
    /** @brief M. Biot alpha constant */
    REAL BiotAlpha()
    {
        return fBiotAlpha;
    }
    
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
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief Print Method */
    virtual void Print(std::ostream &out);
    
private:
    
    /**
     * @brief ZDeformation of the mesh
     */
    STATE fZDeformation;
    
    /** @brief M. Biot Parameter */
    REAL fBiotAlpha;
};



#endif /* defined(__PZ__pzelasticSest2D__) */
