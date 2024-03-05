#ifndef _TPZELASTICITY2DOTITOPO_H_
#define _TPZELASTICITY2DOTITOPO_H_

#include <Elasticity/TPZElasticity2D.h>
#include <TPZMatWithMem.h>


/*
 desired memory:
 sol = solution at the integration point
 dsol = solution's derivatives
 */

//! Data to be stored at each integration point
struct TPZOtiTopoDensity{
    
    STATE fDen = -1.;

    [[nodiscard]] int ClassId() const;
    //! Write to stream(serialization method)
    void Write(TPZStream &buf, int withclassid) const;
    //! Read from stream(serialization method)
    void Read(TPZStream &buf, void *context);
    //! Print contents(debugging method)
    void Print(std::ostream &out) const;
};


inline std::ostream& operator<<( std::ostream& out, const TPZOtiTopoDensity& t ){
    t.Print(out);
    return out;
}

//! Implements a source for scattering problems of planar waveguides
class TPZElasticity2DOtiTopo : public TPZElasticity2D,
public TPZMatWithMem<TPZOtiTopoDensity>{
public:
    
    TPZElasticity2DOtiTopo(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress = 1);
    
    //! All constructors from base class shall be available
//    using TPZElasticity2DOtiTopo::TPZElasticity2DOtiTopo;
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
    TPZElasticity2DOtiTopo * NewMaterial() const override;
    //! Returns name of the class
    std::string Name() const override { return "TPZElasticity2DOtiTopo"; }
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


#endif /* _TPZELASTICITY2DOTITOPO_H_ */
