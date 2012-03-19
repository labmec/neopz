/**
 * @file
 * @brief Contains the TPZFStructMatrix class which implements Full Structural Matrices.
 */

#ifndef TPZFSTRUCTMATRIX_H
#define TPZFSTRUCTMATRIX_H

#include "pzstrmatrix.h"

class TPZCompMesh;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

/**
 * @brief Implements Full Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZFStructMatrix : public TPZStructMatrix {
public:    
	
    TPZFStructMatrix(TPZCompMesh *);
	
    TPZFStructMatrix(TPZAutoPointer<TPZCompMesh> );
    
    virtual TPZMatrix<REAL> * Create();
	
    virtual TPZMatrix<REAL> * CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
public:
	
};

#endif //TPZFSTRUCTMATRIX_H
