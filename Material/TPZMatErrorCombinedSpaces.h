/**
    \file TPZMatErrorCombinedSpaces.h
    Defines the error interface for materials combining multiple approximation spaces.
*/

#ifndef TPZMATERRORCOMBINEDSPACES_H
#define TPZMATERRORCOMBINEDSPACES_H
#include "TPZMatError.h"

template<class T>
class TPZMaterialDataT;
template<class T>
class TPZVec;


template<class TVar>
class TPZMatErrorCombinedSpacesBC;

/**
 * @ingroup combinedspacesinterface
 * @brief Interface for error computation of materials using combined approximation spaces.
 * This material allows the computation of error measures of the FEM approximation based on an exact solution. 
 * It should be used as a template parameter of TPZMatBase.
 */
template<class TVar>
class TPZMatErrorCombinedSpaces : public virtual TPZMatError<TVar>{
 public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatErrorCombinedSpacesBC<TVar>;
    //! Default constructor
    TPZMatErrorCombinedSpaces() = default;

    
//    [[nodiscard]] virtual int ClassId() const override;

    //! @name Error
    /** @{*/
    /*!
      \brief Calculates the error at a given point x.
      \param[in] datavec input data
      \param[out] errors The calculated errors.
     */
    virtual void Errors(const TPZVec<TPZMaterialDataT<TVar>> &data,
                       TPZVec<REAL> &errors) = 0;

    /** @brief Returns the number of error norms.
        Default is 3: norm, L2 and seminorm. */
//    int NEvalErrors() const override {return TPZMatError<TVar>::NEvalErrors();}
    /** @}*/
};

class TPZMaterial;

template<class TVar>
class TPZMatErrorCombinedSpacesBC : public TPZMatErrorCombinedSpaces<TVar>{
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat) {}
public:
    void Errors(const TPZVec<TPZMaterialDataT<TVar>> &data,
                TPZVec<REAL> &errors) override{}
};

extern template class TPZMatErrorCombinedSpaces<STATE>;
extern template class TPZMatErrorCombinedSpaces<CSTATE>;
extern template class TPZMatErrorCombinedSpacesBC<STATE>;
extern template class TPZMatErrorCombinedSpacesBC<CSTATE>;

#endif
