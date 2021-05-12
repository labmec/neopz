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
    void SetMaterialImpl(TPZMaterial* mat){}
};
#endif