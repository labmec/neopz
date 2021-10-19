/**
    \file TPZMatBase.h
    Defines TPZMatBase, which  is the most basic class that 
    materials can inherit from.
*/

#ifndef TPZMATERIALT_H
#define TPZMATERIALT_H

#include "TPZMaterial.h"

template<class TVar>
class TPZMaterialDataT;
template<class TVar>
class TPZBndCondT;


/**
 * @ingroup material
 * @brief This class provides additional interfaces for materials.
 * While the TPZMaterial class had only the type-agnostic interface,
 * this class has type information.
 * @note Actual materials should derive from TPZMatBase
 */
template<class TVar>
class TPZMaterialT : public TPZMaterial{
 public:
    //! Default constructor
    TPZMaterialT() = default;
    //! Constructor taking material identifier
    explicit TPZMaterialT(int id) : TPZMaterial(id){
    }

    /**
       @name ForcingFunction
    */
    /**@{*/
     /**
     * @brief Sets a forcing function to the material.
     * This method is used for setting functions used in the weak formulation.
     * @param[in] f Forcing function
     * @param[in] pOrder Suggested integration rule order.
     */
    virtual void SetForcingFunction(ForcingFunctionType<TVar> f, int pOrder);
    //! Whether a forcing function has been set.
    [[nodiscard]] bool HasForcingFunction() const override{ return (bool)fForcingFunction;}
    //! Get forcing function polynomial order
    [[nodiscard]] int ForcingFunctionPOrder() const
    {
        return fForcingFunctionPOrder;
    }
    //! Gets reference to the forcing function.
    [[nodiscard]] virtual const ForcingFunctionType<TVar> &ForcingFunction() const
    {
        return fForcingFunction;
    }
    [[nodiscard]] virtual ForcingFunctionType<TVar> &ForcingFunction()
    {
        return fForcingFunction;
    }
    /**@}*/
    
    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    virtual TPZBndCondT<TVar>* CreateBC(TPZMaterial *reference,
                                        int id, int type,
                                        const TPZFMatrix<TVar> &val1,
                                        const TPZVec<TVar> &val2);
    
    [[nodiscard]] int ClassId() const override;
 protected:
    ForcingFunctionType<TVar> fForcingFunction{nullptr};
    int fForcingFunctionPOrder{0};
};

extern template class TPZMaterialT<STATE>;
extern template class TPZMaterialT<CSTATE>;
#endif
