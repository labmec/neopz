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

struct TPZElementMatrix;
class TPZCompMesh;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSpStructMatrix : public TPZStructMatrix {
public:    
	
    TPZSpStructMatrix(TPZCompMesh *);
	
    virtual TPZMatrix<REAL> * Create();
	
	virtual TPZMatrix<REAL> * CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone(); 	
	
    /** Used only for testing */
	static int main();
	
};

#endif //TPZSPSTRUCTMATRIX_H
