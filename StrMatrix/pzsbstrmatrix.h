/**
 * @file
 * @brief Contains the TPZSBandStructMatrix class which implements Symmetric Banded Structural Matrices.
 */

#ifndef TPZSBANDSTRUCTMATRIX_H
#define TPZSBANDSTRUCTMATRIX_H

#include "pzstrmatrix.h"

class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;

/**
 * @brief Implements Symmetric Banded Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSBandStructMatrix : public TPZStructMatrix {
public:    
	
    TPZSBandStructMatrix(TPZCompMesh *);
	
    virtual TPZMatrix * Create();
	
    virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
public:
	
};

#endif //TPZSBANDSTRUCTMATRIX_H
