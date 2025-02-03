#ifndef _TPZELASTICITY3DOTITOPO_H_
#define _TPZELASTICITY3DOTITOPO_H_

#include <Elasticity/TPZElasticity3D.h>
#include "TPZElasticity2DOtiTopo.h"
#include <TPZMatWithMem.h>

//! Implements a source for scattering problems of planar waveguides
class TPZElasticity3DOtiTopo : public TPZElasticity3D,
public TPZMatWithMem<TPZOtiTopoDensity>{
public:
    
    TPZElasticity3DOtiTopo(int id, STATE E, STATE nu, TPZVec<STATE> &force);
    
    //! All constructors from base class shall be available
//    using TPZElasticity3DOtiTopo::TPZElasticity3DOtiTopo;
    //! Contribution to the matrix and rhs at the integration point
//    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
//                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override
//    {
//        Contribute(data,weight,ef);
//    }
//    //! Contribution to the rhs at the integration point
//    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
//                    TPZFMatrix<STATE> &ef) override;
    
//    void FillDataRequirements(TPZMaterialData &data) const override;

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void Solution(const TPZMaterialDataT<STATE> &data, int var,
                  TPZVec<STATE> &Solout) override;
    
    //! Creates a copy of this instance
    TPZElasticity3DOtiTopo * NewMaterial() const override;
    //! Returns name of the class
    std::string Name() const override { return "TPZElasticity3DOtiTopo"; }
    //! Unique identifier for serialization purposes
    [[nodiscard]] int ClassId() const override;
    //! Write to stream(serialization method)
    void Write(TPZStream &buf, int withclassid) const override;
    //! Read from stream(serialization method)
    void Read(TPZStream &buf, void *context) override;
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
    
    /** @brief Calculates the element stiffness matrix */
    void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight,
                    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    
    /**
     * @brief Returns the number of variables associated with the variable indexed by var.
     */
    int NSolutionVariables(int var) const override;
};


#endif /* _TPZElasticity3DOtiTopo_H_ */
