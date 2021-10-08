#ifndef TPZMATERIALDATAT_H
#define TPZMATERIALDATAT_H
#include "TPZMaterialData.h"

//! Alias for the type of FEM solution
template<class TVar>
using TPZFemSol = TPZManVector<TVar, TPZMaterialData::MatDataDimSol>;

//! Alias for the type of a FEM solution gradient
template<class TVar>
using TPZFemGradSol = TPZFNMatrix<3*TPZMaterialData::MatDataDimSol, TVar>;
//! Alias for the type of a vector of FEM solutions
template<class TVar>
using TPZSolVec = TPZManVector<TPZFemSol<TVar>,TPZMaterialData::MatDataNumSol>;
//!Alias for the type of a vector of FEM solution gradients
template<class TVar>
using TPZGradSolVec = TPZManVector<TPZFemGradSol<TVar>,TPZMaterialData::MatDataNumSol>;

/**
 * @ingroup material
 * @brief This class implements an interface between TPZCompEl::CalcStiff and Contribute methods of the materials.
 * Since it takes the solution type into account, it can store the solution at an integration point.
 */
template<class TVar>
class TPZMaterialDataT : public TPZMaterialData{
public:
    //! Default constructor.
    TPZMaterialDataT();
    /** @name Data
     * @brief Attributes to be computed in CalcStiff.
     * @{ */
    //! Vector of the solutions at the integration point
    TPZSolVec<TVar> sol;
    //! Vector of the derivatives of the solution at the integration point
    TPZGradSolVec<TVar> dsol;
    //! Vector of the divergence of the solution at the integration point (only of H(div) spaces)
    TPZSolVec<TVar> divsol;
    //! Vector of the curl of the solution at the integration point (only of H(curl) spaces).
    TPZSolVec<TVar> curlsol;
    /** @brief Computes the flux divergence values based on a Material of H(div) approx space */
    /**@}*/
    [[deprecated("This call is useless")]] void ComputeFunctionDivergence() override;

    void SetSolSizes(const int nSol, const int uLen,
                     const int duRow, const int duCol) override;
protected:
    //! DummyFunction to force the parent class to be abstract
    bool HasSol() override{ return true;}
};

extern template class TPZMaterialDataT<STATE>;
extern template class TPZMaterialDataT<CSTATE>;
#endif