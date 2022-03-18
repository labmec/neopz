/**
    \file TPZMatBase.h
    Defines TPZMatBase, which  is the most basic class that 
    materials can inherit from.
*/

#ifndef TPZMATBASE_H
#define TPZMATBASE_H

#include "TPZMaterialT.h"
#include "TPZBndCondBase.h"
#include "pzvec.h"
#include "Hash/TPZHash.h"
/*!
  \brief The class that Materials should inherit from.
  One mandatory interface is either TPZMatSingleSpaceT, for formulations with only one 
  approximation space, or TPZMatCombinedSpacesT, for formulations that use multiple 
  approximation spaces. The other interfaces will enhance the functionality of the material.
  \ingroup material
 */
template<class TVar, class... Interfaces>
class TPZMatBase : public TPZMaterialT<TVar>,
                   public Interfaces... {
 public:
    //! Default constructor
    TPZMatBase() = default;
    //! Constructor taking material identifier
    explicit TPZMatBase(int id) : TPZMaterialT<TVar>(id){
    }

    /** @brief Creates an object TPZBndCond derived of TPZMaterial.
     * The TPZBndCond is a material for calculating boundary conditions.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    TPZBndCondT<TVar>* CreateBC(TPZMaterial *reference,
                         int id, int type,
                         const TPZFMatrix<TVar> &val1,
                         const TPZVec<TVar> &val2) override;

    template<class... AddInterfaces>
    static TPZBndCondT<TVar> *CreateBCImpl(TPZMaterial *reference,
                                    int id, int type,
                                    const TPZFMatrix<TVar> &val1,
                                    const TPZVec<TVar> &val2);


    [[nodiscard]] int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

    
    /** @name MandatoryOverride*/
    /**@{*/
    [[nodiscard]] int Dimension() const override = 0;
    [[nodiscard]] int NStateVariables() const override= 0;
    /**@}*/
    /** @name OptionalOverride*/
    /**@{*/
    /** @brief Returns the name of the material.*/
    [[nodiscard]] std::string Name() const override
    {
        return TPZMaterialT<TVar>::Name();
    }
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    [[nodiscard]] int NSolutionVariables(int var) const override
    {
        return TPZMaterialT<TVar>::NSolutionVariables(var);
    }
    /** @brief Returns the variable index associated with a given name */
    [[nodiscard]] int VariableIndex(const std::string &name) const override
    {
        return TPZMaterialT<TVar>::VariableIndex(name);
    }
    /**@}*/
 protected:
    bool IsMatImpl() final{return true;}
};

// http://labmec.github.io/neopz/material/index.html#boundary-conditions

template<class TVar, class... Interfaces>
template<class... AddInterfaces>
TPZBndCondT<TVar> *
TPZMatBase<TVar, Interfaces...>::CreateBCImpl(TPZMaterial *reference, int id, int type, const TPZFMatrix<TVar> &val1,
                                              const TPZVec<TVar> &val2) {
    return new TPZBndCondBase<TVar, typename Interfaces::TInterfaceBC...,
            typename AddInterfaces::TInterfaceBC...>(reference, id, type, val1, val2);
}

template<class TVar, class... Interfaces>
TPZBndCondT<TVar> *TPZMatBase<TVar, Interfaces...>::CreateBC(TPZMaterial *ref, int id, int type,
                                                             const TPZFMatrix<TVar> &val1, const TPZVec<TVar> &val2) {
    return CreateBCImpl<>(ref, id, type, val1, val2);
}

template<class TVar, class...Interfaces>
int TPZMatBase<TVar, Interfaces...>::ClassId() const{
    constexpr int nInterfaces = sizeof...(Interfaces);
    int id = Hash("TPZMatBase") ^ ClassIdOrHash<TVar>() ^
        TPZMaterial::ClassId() << 1;
    
//    if constexpr (nInterfaces>0){
//        TPZVec<int> classIds = {Interfaces::ClassId()...};
//        for(auto count =0;count< nInterfaces;count++){
//            id = id ^ (classIds[count] << (count+2));
//        }
//    }
    return id;
}

template<class TVar, class...Interfaces>
void TPZMatBase<TVar, Interfaces...>::Read(
    TPZStream& buf, void* context)
{
    TPZMaterialT<TVar>::Read(buf,context);
    /* The following will perform calls to all the Read methods of the interfaces. 
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Read(buf,context),...);
}

template<class TVar, class...Interfaces>
void TPZMatBase<TVar, Interfaces...>::Write(
    TPZStream& buf, int withclassid) const
{
    TPZMaterialT<TVar>::Write(buf,withclassid);
    /* The following will perform calls to all the Read methods of the interfaces. 
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Write(buf,withclassid),...);
}

#endif
