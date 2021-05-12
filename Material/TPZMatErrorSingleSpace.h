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
class TPZMatErrorSingleSpace : public TPZMatError<TVar>{
 public:
    using TInterfaceBC = TPZMatErrorSingleSpaceBC<TVar>;
    //! Default constructor
    TPZMatErrorSingleSpace() = default;

    [[nodiscard]] int ClassId() const override;
    //!@name Error
    /** @{*/
    //! Public interface for calculating errors
    void Errors(const TPZMaterialDataT<TVar> &data,
                TPZVec<REAL> &errors);
    /** @brief Returns the number of error norms.
        Default is 3: norm, L2 and seminorm. */
    int NEvalErrors() override{return TPZMatError<TVar>::NEvalErrors();}
    /** @}*/
protected:
    /*!
      \brief Calculates the error at a given point x.
      \param[in] x coordinates in the domain (not reference element)
      \param[in] sol Solution at point `x`
      \param[in] dsol Solution derivative  at point `x`
      \param[in] axes Transformation between element axis and R3 space
      \param[in] errors The calculated errors.
     */
    virtual void Errors(const TPZVec<REAL> &x, const TPZVec<TVar> &sol,
                        const TPZFMatrix<TVar> &dsol,
                        const TPZFMatrix<REAL> &axes,
                        TPZVec<REAL> &errors) = 0;
};

template<class TVar>
class TPZMatErrorSingleSpaceBC : public TPZMatErrorSingleSpace<TVar>{
protected:
    void SetMaterialImpl(TPZMaterial *mat) {}
    void Errors(const TPZVec<REAL> &x, const TPZVec<TVar> &sol,
                const TPZFMatrix<TVar> &dsol,
                const TPZFMatrix<REAL> &axes,
                TPZVec<REAL> &errors) override{};
};
extern template class TPZMatErrorSingleSpace<STATE>;
extern template class TPZMatErrorSingleSpace<CSTATE>;
extern template class TPZMatErrorSingleSpaceBC<STATE>;
extern template class TPZMatErrorSingleSpaceBC<CSTATE>;
#endif