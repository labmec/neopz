
#include "pzskylstrmatrix.h"
#include "pzcmesh.h" 
#include "pzelmat.h"

#ifndef TPZPARSKYLINESTRUCTMATRIX_H
#define TPZPARSKYLINESTRUCTMATRIX_H


class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;
/**
     Defines Parallel Structural Matrix
     @ingroup structural
*/
class TPZParSkylineStructMatrix : public TPZSkylineStructMatrix {
public:    
    static int main();
    
    TPZParSkylineStructMatrix(TPZCompMesh *);
    
    TPZParSkylineStructMatrix(const TPZParSkylineStructMatrix &cp);

    virtual TPZMatrix * Create();

    virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs);

    virtual TPZStructMatrix * Clone();

public:
};


#endif

