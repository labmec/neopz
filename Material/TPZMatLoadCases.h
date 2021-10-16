#ifndef TPZMATLOADCASES_H
#define TPZMATLOADCASES_H

#include "TPZSavable.h"

class TPZStream;
/*!
  \brief Defines a type-agnostic interface for materials with multiple load cases.
  This class is only useful for dynamic casting.
  @note The template parameter of TPZMatLoadCases is relevant on TPZMatLoadCasesBC
  \ingroup material
 */
class TPZMatLoadCasesBase : public virtual TPZSavable{
public:
    //! Default constructor
    TPZMatLoadCasesBase() = default;
    //!@name LoadCases
    /** @{*/
    //! Returns the number of load cases (rhs columns).
    [[nodiscard]] int NumLoadCases() const {return fNumLoadCases;}
    //! Returns the minimum number of load cases (rhs columns).
    [[nodiscard]] virtual int MinimumNumberofLoadCases() const { return 1;}
    //! Sets the number of load cases (rhs columns).
    void SetNumLoadCases(int numloadcases);
    /** @brief Indicates which variable should be post processed.
     * The method Solution in the material should then pick the 
     * corresponding solution accordingly. 
     */
    void SetPostProcessIndex(int index);
    /** @}*/

    //!@name ReadWrite
    /** @{*/
    //! Read and Write methods
    [[nodiscard]] int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    /** @}*/
protected:
    //! Number of expected columns in the rhs vector.
    int fNumLoadCases{1};
    //! Indicates with solution should be used for post-processing.
    int fPostProcIndex{0};
    //! Dummy method for preventing its creation
    virtual void AmIConcrete() const = 0;
};

//forward declaration
template<class TVar>
class TPZMatLoadCasesBC;
/*!
  \brief Defines an interface for materials with multiple load cases.
  Materials with this interface are able to solve problems in which the
  right hand side vector (and thus the solution) has more than one column.
  @note The template parameter  is relevant on TPZMatLoadCasesBC
  \ingroup singlespaceinterface combinedspacesinterface
 */
template<class TVar>
class TPZMatLoadCases : public TPZMatLoadCasesBase{
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatLoadCasesBC<TVar>;
protected:
    //! Dummy method for ensuring that its base class is not used
    void AmIConcrete() const final {}
};

#include "pzvec.h"

class TPZMaterial;
template<class T>
class TPZFMatrix;

/**
   @brief Company class for TPZMatLoadCases.
   This class allows boundary conditions to define several values for the right hand side value.
 */
template<class TVar>
class TPZMatLoadCasesBC : public TPZMatLoadCases<TVar>{
protected:
    TPZMatLoadCases<TVar> *fMatLoadCases{nullptr};
    TPZVec<TPZVec<TVar>> fBCRhsValVec;
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat);
public:
    //! Set a vector boundary condition rhs values
    void SetBCRhsValVec(TPZVec<TPZVec<TVar>>& bcValVec);
    //! Get the `i`th boundary condition rhs value
    const TPZVec<TVar>& GetBCRhsVal(int i) const;
};

extern template class TPZMatLoadCases<STATE>;
extern template class TPZMatLoadCases<CSTATE>;
extern template class TPZMatLoadCasesBC<STATE>;
extern template class TPZMatLoadCasesBC<CSTATE>;
#endif
