/**
    \file TPZMatErrorSingleSpace.h
    Defines the error interface for single approximation space materials.
*/

#ifndef TPZMATERRORSINGLESPACE_H
#define TPZMATERRORSINGLESPACE_H
#include "TPZMatError.h"

#include "TPZMaterialDataT.h"


template<class TVar>
class TPZMatErrorSingleSpaceBC;
/**
 * @ingroup singlespaceinterface
 * @brief Interface for error computation of materials using a single approximation space.
 * This material allows the computation of error measures of the FEM approximation based on an exact solution. 
 */
template<class TVar>
class TPZMatErrorSingleSpace : public virtual TPZMatError<TVar>{
 public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatErrorSingleSpaceBC<TVar>;
    //! Default constructor
    TPZMatErrorSingleSpace() = default;

//    [[nodiscard]] int ClassId() const override;
    //!@name Error
    /** @{*/
    /** @brief Public interface for calculating errors
     See TPZMaterialData and TPZMaterialDataT for the available data.*/
    virtual void Errors(const TPZMaterialDataT<TVar> &data,
                TPZVec<REAL> &errors) = 0;
    /** @brief Returns the number of error norms.
        Default is 3: norm, L2 and seminorm. */
//    int NEvalErrors() const override{return TPZMatError<TVar>::NEvalErrors();}
    /** @}*/
};

//!BC Interface for TPZMatErrorSingleSpace
template<class TVar>
class TPZMatErrorSingleSpaceBC : public TPZMatErrorSingleSpace<TVar>{
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat) {}
    void Errors(const TPZMaterialDataT<TVar> &data,
                TPZVec<REAL> &errors) override{};
};
extern template class TPZMatErrorSingleSpace<STATE>;
extern template class TPZMatErrorSingleSpace<CSTATE>;
extern template class TPZMatErrorSingleSpaceBC<STATE>;
extern template class TPZMatErrorSingleSpaceBC<CSTATE>;
#endif
