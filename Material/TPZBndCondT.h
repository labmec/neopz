#ifndef TPZBNDCONDT_H
#define TPZBNDCONDT_H

#include "TPZBndCond.h"
#include "TPZMatTypes.h"
#include "pzfmatrix.h"

template<class T>
class TPZMaterialDataT;
/**
 * @ingroup matsinglespace
 * @brief This class defines an interface for boundary conditions.
 * It complements the interfaces defined in TPZBndCond with type information.
 * @note Boundary conditions will actually be created as TPZBndCondBase
 */
template <class TVar>
class TPZBndCondT : public TPZBndCond
{
protected:
    //! Default constructor
    TPZBndCondT() = default;
    //! Constructor setting type, fVal1 and fVal2
    TPZBndCondT(int type,
                const TPZFMatrix<TVar> &val1,
                const TPZVec<TVar> &val2);
public:
    /*no need to specify destructor, copy/move assignment/destructors.
     the compiler can do it for us*/

    /**
     * @brief Sets a forcing function to the boundary material.
     * This method is used for setting functions used in the weak formulation.
     * @param[in] f Forcing function
     * @param[in] pOrder Suggested integration rule order.
     */
    void SetForcingFunctionBC(ForcingFunctionBCType<TVar> f){
        fForcingFunctionBC = f;
    }
    //! Whether a forcing function has been set.
    [[nodiscard]] bool HasForcingFunctionBC() const final{
        return (bool)fForcingFunctionBC;}
    //@{
    //! Gets reference to the boundary condition forcing function.
    [[nodiscard]] const ForcingFunctionBCType<TVar> &ForcingFunctionBC() const
    {
        return fForcingFunctionBC;
    }
    [[nodiscard]] ForcingFunctionBCType<TVar> &ForcingFunctionBC()
    {
        return fForcingFunctionBC;
    }
    //@}

    /** @name ReadWrite */
	/** @{ */
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    /** @} */
    //! Gets the FEM matrix value for the boundary condition.
    const TPZFMatrix<TVar> &Val1() const
    { return fBCVal1; }
	//! Gets the right hand side value for the boundary condition.
	const TPZVec<TVar> &Val2() const
    { return fBCVal2; }
    //! Sets the FEM matrix value for the boundary condition.
    void SetVal1(const TPZFMatrix<TVar> &val1)
    { fBCVal1 = val1; }
	//! Sets the right hand side value for the boundary condition.
	void SetVal2(const TPZVec<TVar> &val2)
    { fBCVal2 = val2; }
protected:
    /** @brief first value (FEM matrix)  of boundary condition */
 	TPZFNMatrix<9,TVar> fBCVal1;
	/** @brief second value (rhs vector) of boundary condition */
	TPZManVector<TVar,3> fBCVal2;
    /** @brief Boundary condition forcing function */
    ForcingFunctionBCType<TVar> fForcingFunctionBC;
};

extern template class TPZBndCondT<STATE>;
extern template class TPZBndCondT<CSTATE>;
#endif