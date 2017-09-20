/**
 * @file
 * @brief Contains the TPZParSkylineStructMatrix class which defines parallel structural matrix for skyline matrices.
 */

#ifndef TPZPARSKYLINESTRUCTMATRIX_H
#define TPZPARSKYLINESTRUCTMATRIX_H

#include "pzskylstrmatrix.h"
#include "pzcmesh.h" 
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"


/**
 * @brief Defines parallel structural matrix for skyline matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZParSkylineStructMatrix : public TPZSkylineStructMatrix {
	
	int fNumThreads;
	
public:
	static int main();
	
	TPZParSkylineStructMatrix(TPZCompMesh *, int numthreads);
    
    TPZParSkylineStructMatrix(const TPZParSkylineStructMatrix &cp);
	
    virtual TPZMatrix<STATE> * Create();
	
    using TPZStructMatrix::CreateAssemble;
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
public:
	
};

#endif

