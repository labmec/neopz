/**
 * @file
 * @brief Contains the TPZSymetricSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSymetricSpStructMatrix_H
#define TPZSymetricSpStructMatrix_H

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
class TPZSymetricSpStructMatrix : public TPZStructMatrix {
    
public:    
	
    TPZSymetricSpStructMatrix(TPZCompMesh *);
	
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZMatrix<STATE> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    
    using TPZStructMatrix::CreateAssemble;
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
    
	
    /** Used only for testing */
	static int main();
	
private :
    
    TPZSymetricSpStructMatrix();
    
    friend TPZPersistenceManager;
};

#endif //TPZSymetricSpStructMatrix_H
