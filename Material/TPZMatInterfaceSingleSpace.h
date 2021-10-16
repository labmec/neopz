//
// Created by Gustavo Batistela on 4/28/21.
//

#ifndef TPZMATINTERFACESINGLESPACE_H
#define TPZMATINTERFACESINGLESPACE_H

#include "TPZSavable.h"
#include "TPZMaterialDataT.h"
template<class TVar>
class TPZBndCondT;
template<class Tvar>
class TPZMatInterfaceSingleSpaceBC;


/**
 * @ingroup singlespaceinterface
 * @brief Interface for discontinuous materials using a single approximation space.
 * This material allows contribution of the FEM matrix using the interface between elements. It should be used as a template parameter of TPZMatBase.
 */
template<class TVar>
class TPZMatInterfaceSingleSpace : public virtual TPZSavable {
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatInterfaceSingleSpaceBC<TVar>;
    [[nodiscard]] int ClassId() const override;
    //!@name Interface
    /** @{*/
    /**
     * @brief Returns the solution associated with the var index based on the finite element approximation
     */
    virtual void
    SolutionInterface(const TPZMaterialDataT<TVar> &data,
                      const TPZMaterialDataT<TVar> &dataleft,
                      const TPZMaterialDataT<TVar> &dataright,
                      const int var,
                      TPZVec<TVar> &Solout) = 0;

    /**
     * @brief It computes a contribution to stiffness matrix and load vector at one integration point
     * @param [in] data
     * @param [in] dataleft
     * @param [in] dataright
     * @param [in] weight
     * @param [out] ek is the stiffness matrix
     * @param [out] ef is the load vector
     * @since April 16, 2007
     */
    virtual void
    ContributeInterface(const TPZMaterialDataT<TVar> &data,
                        const TPZMaterialDataT<TVar> &dataleft,
                        const TPZMaterialDataT<TVar> &dataright,
                        REAL weight, TPZFMatrix<TVar> &ek,
                        TPZFMatrix<TVar> &ef) = 0;

    /**
     * @brief It computes a contribution to residual vector at one integration point
     * @param [in] data
     * @param [in] dataleft
     * @param [in] dataright
     * @param [in] weight
     * @param [out] ef is the load vector
     * @since April 16, 2007
     */
    virtual void
    ContributeInterface(const TPZMaterialDataT<TVar> &data,
                        const TPZMaterialDataT<TVar> &dataleft,
                        const TPZMaterialDataT<TVar> &dataright,
                        REAL weight, TPZFMatrix<TVar> &ef);

    /**
     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
     * @param [in] data
     * @param [in] dataleft
     * @param [in] weight
     * @param [out] ek is the stiffness matrix
     * @param [out] ef is the load vector
     * @param [in] bc is the boundary condition object
     * @since February 21, 2013
     */
    virtual void
    ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                          const TPZMaterialDataT<TVar> &dataleft, REAL weight,
                          TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                          TPZBndCondT<TVar> &bc) = 0;


    /**
     * @brief It computes a contribution to residual vector at one BC integration point
     * @param [in] data
     * @param [in] dataleft
     * @param [in] weight
     * @param [out] ef is the load vector
     * @param [in] bc is the boundary condition object
     * @since April 16, 2007
     */
    virtual void
    ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                          const TPZMaterialDataT<TVar> &dataleft, REAL weight,
                          TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc);

    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution
     * of interface elements
     */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data) const = 0;
    /** @}*/
};

template<class TVar>
class TPZMatInterfaceSingleSpaceBC :
    public TPZMatInterfaceSingleSpace<TVar>{
protected:
    TPZMatInterfaceSingleSpace<TVar>* fMatInterface{nullptr};
        // this method is your chance to verify if the material to which this
        // BC interface applies is compatible with this boundary interface
        // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat);
public:    
    void
    SolutionInterface(const TPZMaterialDataT<TVar> &data,
                      const TPZMaterialDataT<TVar> &dataleft,
                      const TPZMaterialDataT<TVar> &dataright,
                      const int var,
                      TPZVec<TVar> &Solout) override;

     void
    ContributeInterface(const TPZMaterialDataT<TVar> &data,
                        const TPZMaterialDataT<TVar> &dataleft,
                        const TPZMaterialDataT<TVar> &dataright,
                        REAL weight, TPZFMatrix<TVar> &ek,
                        TPZFMatrix<TVar> &ef) override;

    void
    ContributeInterface(const TPZMaterialDataT<TVar> &data,
                        const TPZMaterialDataT<TVar> &dataleft,
                        const TPZMaterialDataT<TVar> &dataright,
                        REAL weight, TPZFMatrix<TVar> &ef) override;
    //throws
    void
    ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                          const TPZMaterialDataT<TVar> &dataleft, REAL weight,
                          TPZFMatrix<TVar> &ef,
                          TPZBndCondT<TVar> &bc) override;
    //throws
    void
    ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                          const TPZMaterialDataT<TVar> &dataleft, REAL weight,
                          TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                          TPZBndCondT<TVar> &bc) override;
    
    void FillDataRequirementsInterface(TPZMaterialData &data) const override;
};

extern template class TPZMatInterfaceSingleSpace<STATE>;
extern template class TPZMatInterfaceSingleSpace<CSTATE>;
extern template class TPZMatInterfaceSingleSpaceBC<STATE>;
extern template class TPZMatInterfaceSingleSpaceBC<CSTATE>;
#endif //TPZMATINTERFACESINGLESPACE_H
