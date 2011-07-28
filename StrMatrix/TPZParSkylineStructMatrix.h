
#include "pzskylstrmatrix.h"
#include "pzcmesh.h" 
#include "pzelmat.h"

#ifndef TPZPARSKYLINESTRUCTMATRIX_H
#define TPZPARSKYLINESTRUCTMATRIX_H


class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;
/**
     @brief Defines Parallel Structural Matrix
     @ingroup structural
*/
class TPZParSkylineStructMatrix : public TPZSkylineStructMatrix {

	int fNumThreads;

public:
		static int main();

		TPZParSkylineStructMatrix(TPZCompMesh *, int numthreads);
    
    TPZParSkylineStructMatrix(const TPZParSkylineStructMatrix &cp);

    virtual TPZMatrix * Create();

    virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    virtual TPZStructMatrix * Clone();

public:
};


#endif

