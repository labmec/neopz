/**
 * @file
 * @brief Contains the TPZElasticity3D class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef TPZHYBRIDELASTICITY3D_H
#define TPZHYBRIDELASTICITY3D_H


#include "TPZElasticity3D.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZHybridElasticity3D : public TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>, public TPZElasticity3D
{
    using TBase = TPZMatBase<STATE,
                             TPZMatSingleSpaceT<STATE>,
                             TPZMatErrorSingleSpace<STATE>,
                             TPZMatLoadCases<STATE>>;
public :
    /**
     * @brief Creates an elastic material with:
     * @param id material id
     * @param E elasticity modulus
     * @param nu poisson coefficient
     * @param fx forcing function \f$ -x = fx \f$
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    TPZHybridElasticity3D(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                          STATE preStressXX = 0., STATE preStressYY = 0., STATE preStressZZ = 0.) :
    TPZElasticity3D(nummat, E, poisson, force,
                    preStressXX, preStressYY, preStressZZ)
    {
        
    }
    
    TPZHybridElasticity3D(int id) : TPZElasticity3D(id) {}
    
    TPZHybridElasticity3D() : TPZElasticity3D(){}

    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    virtual TPZBndCondT<STATE>* CreateBC(TPZMaterial *reference,
                                        int id, int type,
                                        const TPZFMatrix<STATE> &val1,
                                        const TPZVec<STATE> &val2) override
    {
        return new  TPZBndCondBase<STATE,TPZMatCombinedSpacesBC<STATE>, TPZMatErrorCombinedSpacesBC<STATE> >
        (reference,id, type,val1,val2);
    }

    /** @brief Returns the material name*/
    std::string Name()  const override { return "TPZHybridElasticity3D"; }
        
    /** @name Contribute methods */
    /** @{ */
    
    /** @brief Calculates the element stiffness matrix */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data, STATE weight,
                    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Applies the element boundary conditions */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,STATE weight,
                      TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /*
     * @brief Fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    
    /** @} */

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &data, int var,
                  TPZVec<STATE> &Solout) override;
    
    /** @} */

    /** @brief Creates a new material from the current object   ??*/
    TPZMaterial * NewMaterial()  const override;
    /**
       @name ReadWrite
       @{*/
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream &buf, void *context) override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    /**@}*/
    /**
     * @brief Computes the error due to the difference between the interpolated flux \n
     * and the flux computed based on the derivative of the solution
     */

    /** @name Errors */
    /** @{ */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                TPZVec<STATE> &values) override;
    /** @} */
    
    
};

#endif

