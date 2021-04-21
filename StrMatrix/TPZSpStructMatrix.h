/**
 * @file
 * @brief Contains the TPZSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSPSTRUCTMATRIX_H
#define TPZSPSTRUCTMATRIX_H

#include "TPZStructMatrix.h"

#include "pzstrmatrixor.h"
/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSpStructMatrix : public TPZStructMatrix,
                                  public TPar{
public:
    TPZSpStructMatrix(TPZCompMesh* m) : TPZStructMatrix(m){}
    TPZSpStructMatrix(TPZAutoPointer<TPZCompMesh> m) :
        TPZStructMatrix(m){}
    
    TPZMatrix<TVar> *Create() override;
	TPZStructMatrix * Clone() override;
    
	TPZMatrix<TVar> *CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
    TPZSpStructMatrix() = default;
    
    friend TPZPersistenceManager;
};

#endif //TPZSPSTRUCTMATRIX_H
