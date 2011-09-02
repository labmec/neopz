/**
 * @file
 * @brief Contains the TPBSpStructMatrix class which assembles on the pair equations.
 */

#ifndef TPBSPSTRUCTMATRIX_H
#define TPBSPSTRUCTMATRIX_H

#include "TPZSpStructMatrix.h"

struct TPZElementMatrix;
class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;
class TPZStructMatrix;

/** 
 * @ingroup structural
 * @brief Assembles on the pair equations. \ref structural "Structural Matrix"
 */
class TPBSpStructMatrix : public TPZSpStructMatrix {
public:    
	
    virtual TPZMatrix * Create();    
	
    virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
    /** Used only for testing */
    static int main();
	
    TPBSpStructMatrix(TPZCompMesh *);
    
    TPBSpStructMatrix(const TPBSpStructMatrix &copy) : TPZSpStructMatrix(copy)
    {
    }
	
};

#endif //TPBSPSTRUCTMATRIX_H
