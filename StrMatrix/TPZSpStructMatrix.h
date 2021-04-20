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
class TPZSpStructMatrix : public TPZStructMatrix,
                          public TPZStructMatrixOR {
public:
    TPZSpStructMatrix(TPZCompMesh* m) : TPZStructMatrix(m){}
    TPZSpStructMatrix(TPZAutoPointer<TPZCompMesh> m) :
        TPZStructMatrix(m){}
    
    TPZBaseMatrix * Create() override;
	TPZStructMatrix * Clone() override;
    
	TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
	
    /** Used only for testing */
	static int main();
private :
    TPZSpStructMatrix() = default;
    
    friend TPZPersistenceManager;
};

#endif //TPZSPSTRUCTMATRIX_H
