/**
 * @file
 * @brief Contains the TPZSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSPSTRUCTMATRIX_H
#define TPZSPSTRUCTMATRIX_H

#include "pzstrmatrix.h"

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSpStructMatrix : public TPZStructMatrix {
public:    
	
    TPZSpStructMatrix(TPZCompMesh *);
	
    virtual TPZBaseMatrix * Create() override;
	
    using TPZStructMatrix::CreateAssemble;
    
	TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
    virtual TPZStructMatrix * Clone() override;
    
    public:
int ClassId() const override;

	
    /** Used only for testing */
	static int main();
private :
    TPZSpStructMatrix();
    
    friend TPZPersistenceManager;
};

#endif //TPZSPSTRUCTMATRIX_H
