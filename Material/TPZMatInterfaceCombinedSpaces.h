//
// Created by Gustavo Batistela on 4/28/21.
//

#ifndef TPZMATINTERFACECOMBINEDSPACES_H
#define TPZMATINTERFACECOMBINEDSPACES_H


#include "TPZSavable.h"
#include "TPZMaterialDataT.h"
template<class TVar>
class TPZBndCondT;

template<class TVar>
class TPZMatInterfaceCombinedSpacesBC;

/**
 * @ingroup combinedspacesinterface
 * @brief Interface for discontinuous materials using combined approximation spaces.
 * This material allows contribution of the FEM matrix using the interface between elements. It should be used as a template parameter of TPZMatBase.
 */
template<class TVar>
class TPZMatInterfaceCombinedSpaces : public virtual TPZSavable {
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatInterfaceCombinedSpacesBC<TVar>;
    [[nodiscard]] int ClassId() const override;
    //!@name Interface
    /** @{*/
    /**
     * @brief returns the integration order as a function of interpolation orders of the left and right elements
     */
    [[nodiscard]] virtual int GetIntegrationOrder(const TPZVec<int> &porder_left, const TPZVec<int> &porder_right) const;

    /**
     * @brief Computes a contribution to the stiffness matrix and load vector at one integration point
     * to multiphysics simulation
     * @param [in] data
     * @param [in] dataleft
     * @param [in] dataright
     * @param [in] weight
     * @param [out] ek is the stiffness matrix
     * @param [out] ef is the load vector
     * @since June 5, 2012
     */
    virtual void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                                     const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                                     const std::map<int, TPZMaterialDataT<TVar>> &dataright, REAL weight,
                                     TPZFMatrix<TVar> &ek,
                                     TPZFMatrix<TVar> &ef) = 0;


    /**
     * @brief Computes a contribution to residual vector at one integration point
     * @param [in] data
     * @param [in] dataleft
     * @param [in] dataright
     * @param [in] weight
     * @param [out] ef is the load vector
     * @since June 5, 2012
     */
    virtual void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                                     const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                                     const std::map<int, TPZMaterialDataT<TVar>> &dataright, REAL weight,
                                     TPZFMatrix<TVar> &ef);

    /**
     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
     * to multiphysics simulation
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
                          const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                          REAL weight,
                          TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc) = 0;

    /**
     * @brief It computes a contribution to the load vector at one BC integration point
     * to multiphysics simulation
     * @param [in] data
     * @param [in] dataleft
     * @param [in] weight
     * @param [out] ef is the load vector
     * @param [in] bc is the boundary condition object
     * @since February 21, 2013
     */
    virtual void
    ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                          const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                          REAL weight,
                          TPZFMatrix<TVar> &ef,
                          TPZBndCondT<TVar> &bc);

    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution
     * of interface elements
     */
    virtual void
    FillDataRequirementsInterface(TPZMaterialDataT<TVar> &data,
                                  std::map<int, TPZMaterialDataT<TVar>> &datavec_left,
                                  std::map<int, TPZMaterialDataT<TVar>> &datavec_right) = 0;

    /**
     * @brief Returns the solution associated with the var index based on the finite element approximation around
     * one interface element
     */
    virtual void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                                   const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                                   const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                                   int var, TPZVec<TVar> &Solout) = 0;

    /**
     * @brief Returns the solution associated with the var index based on the finite element approximation around
     * one interface element
     */
    virtual void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                                   const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                                   const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                                   int var, TPZVec<TVar> &Solout,
                                   TPZCompEl *left,TPZCompEl *right) = 0;
    /** @}*/
};

template<class TVar>
class TPZMatInterfaceCombinedSpacesBC :
    public TPZMatInterfaceCombinedSpaces<TVar>{
protected:
    TPZMatInterfaceCombinedSpaces<TVar>* fMatInterface{nullptr};
        // this method is your chance to verify if the material to which this
        // BC interface applies is compatible with this boundary interface
        // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat);
public:    
    void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                             REAL weight,
                             TPZFMatrix<TVar> &ek,
                             TPZFMatrix<TVar> &ef) override;
    
    void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                             REAL weight,
                             TPZFMatrix<TVar> &ef) override;     
    //throws
    void ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                               const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                               REAL weight,
                               TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                               TPZBndCondT<TVar> &bc) override;
    //throws
    void ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                               const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                               REAL weight,
                               TPZFMatrix<TVar> &ef,
                               TPZBndCondT<TVar> &bc) override;
    
    void
    FillDataRequirementsInterface(TPZMaterialDataT<TVar> &data,
                                  std::map<int, TPZMaterialDataT<TVar>> &datavec_left,
                                  std::map<int, TPZMaterialDataT<TVar>> &datavec_right) override;
    
    void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                           const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                           int var, TPZVec<TVar> &Solout) override;

    
    void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                           const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                           int var, TPZVec<TVar> &Solout,
                           TPZCompEl *left,TPZCompEl *right) override;
};


extern template class TPZMatInterfaceCombinedSpaces<STATE>;
extern template class TPZMatInterfaceCombinedSpaces<CSTATE>;
extern template class TPZMatInterfaceCombinedSpacesBC<STATE>;
extern template class TPZMatInterfaceCombinedSpacesBC<CSTATE>;
#endif //TPZMATINTERFACECOMBINEDSPACES_H
