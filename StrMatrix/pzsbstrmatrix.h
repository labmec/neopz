/**
 * @file
 * @brief Contains the TPZSBandStructMatrix class which implements Symmetric Banded Structural Matrices.
 */

#ifndef TPZSBANDSTRUCTMATRIX_H
#define TPZSBANDSTRUCTMATRIX_H

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzstrmatrix.h"
#include "pzcmesh.h"

/**
 * @brief Implements Symmetric Banded Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSBandStructMatrix : public TPZStructMatrix {
public:    
	
    TPZSBandStructMatrix(TPZCompMesh *);
	
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();

private :
    TPZSBandStructMatrix();
	
    friend TPZPersistenceManager;
};

#endif //TPZSBANDSTRUCTMATRIX_H
