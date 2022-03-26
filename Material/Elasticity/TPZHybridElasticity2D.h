/**
 * @file
 * @brief Contains the TPZElasticity2D class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef TPZHYBRIDELASTICITY2D_H
#define TPZHYBRIDELASTICITY2D_H


#include "TPZElasticity2D.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZHybridElasticity2D : public TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>, public TPZElasticity2D
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
	TPZHybridElasticity2D(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress = 1) :
    TPZElasticity2D(id,E,nu,fx,fy,planestress)
    {
        
    }
    
    TPZHybridElasticity2D(int id) : TPZElasticity2D(id) {}
    
    TPZHybridElasticity2D() : TPZElasticity2D(){}

	/** @brief Returns the material name*/
	std::string Name()  const override { return "TPZHybridElasticity2D"; }
		
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
