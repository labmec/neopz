#ifndef TPZMATQUADRATICEIGENVAL
#define TPZMATQUADRATICEIGENVAL

#include "TPZSavable.h"

class TPZStream;

class TPZMatQuadraticEigenValBC;
/**
  @brief Defines an interface for solving Quadratic Eigenvalue problems.
  * This interface simply adds a method to switch between the matrices K,
  * L and M in the equation -b^2 Mu + i b Lu + Ku = 0,
  * where (b,u) is the eigenpair and i is the imaginary unit.
  * It is suggested that the materials using this interface
  * use a function pointer for switching between the correct contribute method.
 */
class TPZMatQuadraticEigenVal : public virtual TPZSavable{
public:
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatQuadraticEigenValBC;
    //! Default constructor
    TPZMatQuadraticEigenVal() = default;
    //! Enum for setting which matrix is being assembled
    enum EWhichMatrix{
        NDefined=0,///< Not defined
        K,///< K
        L,///< L
        M///< M
    };
    /**
       @name Quadratic
       @{
    */
    //! Sets the current matrix to A
    virtual void SetMatrixK(){ fAssembling=K;}
    //! Sets the current matrix to B
    virtual void SetMatrixL(){ fAssembling=L;}
    //! Sets the current matrix to A
    virtual void SetMatrixM(){ fAssembling=M;}
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

class TPZMatQuadraticEigenValBC : public TPZMatQuadraticEigenVal{
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial* mat){}
};
#endif
