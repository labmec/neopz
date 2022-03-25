//
//  hybridpoissoncollapsed.h
//  PZ
//

#ifndef TPZHybridPoissonCollapsed_hpp
#define TPZHybridPoissonCollapsed_hpp

#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMaterialDataT.h"

/**
 * @ingroup material
 * @author Karolinne Coelho
 * @since 11/17/2020
 * @brief Material to solve a mixed poisson problem 2d by multiphysics simulation with HDivCollapsed spaces
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
class TPZHybridPoissonCollapsed : public TPZMixedDarcyFlow{

// class TPZHybridPoissonCollapsed: virtual public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
//         TPZMatErrorCombinedSpaces<STATE>, TPZIsotropicPermeability>, public TPZMixedDarcyFlow {
    
    REAL fPermeability = 1.;

    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
            TPZMatErrorCombinedSpaces<STATE>, TPZIsotropicPermeability>;

public:
    TPZHybridPoissonCollapsed();

        /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZHybridPoissonCollapsed(int id, int dim);
    
    TPZHybridPoissonCollapsed(const TPZHybridPoissonCollapsed &cp);

    ~TPZHybridPoissonCollapsed() override;

    TPZHybridPoissonCollapsed &operator=(const TPZHybridPoissonCollapsed &copy);

    TPZMaterial *NewMaterial() const override{
        return new TPZHybridPoissonCollapsed(*this);
    }

    void SetIsotropicPermeability(const REAL constant)
    {
        fPermeability = constant;
        this->SetConstantPermeability(constant);
    }

    REAL GetPermeability()
    {
        return fPermeability;
    }
	    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) override;
	
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors);

};

#endif
