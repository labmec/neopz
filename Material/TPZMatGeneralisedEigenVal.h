#ifndef TPZMATGENERALISEDEIGENVAL
#define TPZMATGENERALISEDEIGENVAL

#include "TPZSavable.h"

class TPZStream;

class TPZMatGeneralisedEigenValBC;
/**
  @brief Defines an interface for solving Generalised Eigenvalue problems.
  * This interface simply adds a method to switch between the left (A) and the 
  * right (B) matrices. It is suggested that the materials using this interface
  * use a function pointer for switching between the correct contribute method.
  * See TPZWaveguideModalAnalysis for an example of usage.
 */
class TPZMatGeneralisedEigenVal : public virtual TPZSavable{
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatGeneralisedEigenValBC;
    //! Default constructor
    TPZMatGeneralisedEigenVal() = default;
    //! Enum for setting which matrix is being assembled
    enum EWhichMatrix{
        NDefined=0,///< Not defined
        A,///< A (left)
        B///< B (right)
    };
    /**
       @name Generalised
       @{
    */
    //! Sets the current matrix to A
    virtual void SetMatrixA(){ fAssembling=A;}
    //! Sets the current matrix to B
    virtual void SetMatrixB(){ fAssembling=B;}
    /** @} */
    //!@name ReadWrite
    /** @{*/
    //! Read and Write methods
    [[nodiscard]] int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override{};

    void Write(TPZStream& buf, int withclassid) const override{};
    /** @}*/
protected:
    //! Which matrix is currently being assembled
    EWhichMatrix fAssembling{NDefined};
};

class TPZMaterial;

class TPZMatGeneralisedEigenValBC : public TPZMatGeneralisedEigenVal{
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial* mat){}
};
#endif
