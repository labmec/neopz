/**
 * @file
 * @brief Contains the TPZParSkylineStructMatrix class which defines parallel structural matrix for skyline matrices.
 */

#include "pzskylstrmatrix.h"
#include "pzcmesh.h" 
#include "pzelmat.h"

#ifndef TPZPARSKYLINESTRUCTMATRIX_H
#define TPZPARSKYLINESTRUCTMATRIX_H

class TPZCompMesh;
template<class TVar> 
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

/**
 @brief Defines parallel structural matrix for skyline matrices. \ref structural "Structural Matrix"
 @ingroup structural
 */
class TPZParSkylineStructMatrix : public TPZSkylineStructMatrix {
	
	int fNumThreads;
	
public:
	static int main();
	
	TPZParSkylineStructMatrix(TPZCompMesh *, int numthreads);
    
    TPZParSkylineStructMatrix(const TPZParSkylineStructMatrix &cp);
	
    virtual TPZMatrix<REAL> * Create();
	
    virtual TPZMatrix<REAL> * CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
public:
	
};

#endif

