/**
 * @file
 * @brief Contains the TPZSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSPSTRUCTMATRIX_H
#define TPZSPSTRUCTMATRIX_H

#include "pzstrmatrix.h"
#include "pzysmp.h"

#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSpStructMatrix : public TPZStructMatrix {
public:    
	
    TPZSpStructMatrix(TPZCompMesh *);
	
    virtual TPZMatrix<STATE> * Create();
	
    using TPZStructMatrix::CreateAssemble;
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone(); 	
    
    public:
virtual int ClassId() const;

	
    /** Used only for testing */
	static int main();
private :
    TPZSpStructMatrix();
    
    friend TPZPersistenceManager;
};

#endif //TPZSPSTRUCTMATRIX_H
