/**
 * @file
 * @brief Contains the TPZFStructMatrix class which implements Full Structural Matrices.
 */

#ifndef TPZFSTRUCTMATRIX_H
#define TPZFSTRUCTMATRIX_H

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzstrmatrix.h"

class TPZCompMesh;

/**
 * @brief Implements Full Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZFStructMatrix : public TPZStructMatrix {
public:    
	
    TPZFStructMatrix(TPZCompMesh *);
	
    TPZFStructMatrix(TPZAutoPointer<TPZCompMesh> );
    
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZStructMatrix * Clone();

private :
	
    TPZFStructMatrix();
    
    friend TPZPersistenceManager;
};

#endif //TPZFSTRUCTMATRIX_H
