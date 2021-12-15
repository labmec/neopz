#ifndef TPZMATERROR_H
#define TPZMATERROR_H

#include "TPZSavable.h"
#include "pzreal.h"

#include <functional>

template<class T>
class TPZVec;
template<class T>
class TPZFMatrix;

 //! Alias for exact solution type
template<class TVar>
using ExactSolType = std::function<void (const TPZVec<REAL> &loc, TPZVec<TVar> &result,
                                TPZFMatrix<TVar> &deriv)>;
/*!
  \brief Defines the common error interface used in TPZMatErrorSingleSpace and TPZMatErrorCombinedSpaces
  \ingroup material
 */
template<class TVar>
//class TPZMatError : public virtual TPZSavable{
class TPZMatError {
public:

    //! Default constructor
    TPZMatError() = default;
    //!@name Error
    /** @{*/
    /**
     * @brief Sets a exact solution to the material.
     * This method is used for setting functions used in the weak formulation.
     * @param f Forcing function
     * @param pOrder Suggested integration rule order.
     */
    void SetExactSol(ExactSolType<TVar> f, int pOrder);
    //! Whether a exact solution has been set.
    [[nodiscard]] bool HasExactSol() const { return (bool)fExactSol;}
    //! Gets reference to the exact solution.
    [[nodiscard]] const ExactSolType<TVar> &ExactSol() const
    {
        return fExactSol;
    }
    //! Gets reference to the exact solution.
    [[nodiscard]] ExactSolType<TVar> &ExactSol()
    {
        return fExactSol;
    }
    //! Gets polynomial order of exact solution.
    [[nodiscard]] int PolynomialOrderExact() const
    {return fExactPOrder;}

    /** @brief Returns the number of error norms.
        Default is 3: norm, L2 and seminorm. */
    virtual int NEvalErrors() const{return 3;}
    /** @}*/

//    [[nodiscard]] int ClassId() const;
protected:
    //! Exact Solution
    ExactSolType<TVar> fExactSol = 0;
    //! Suggested integration order for ExactSol
    int fExactPOrder{0};
};

extern template class TPZMatError<STATE>;
extern template class TPZMatError<CSTATE>;
#endif
