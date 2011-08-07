/**
 * @file
 * @brief Contains the TPZSkylineStructMatrix class which implements SkyLine Structural Matrices.
 */

#ifndef TPZSKYLINESTRUCTMATRIX_H
#define TPZSKYLINESTRUCTMATRIX_H

#include "pzstrmatrix.h"

class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;

/**
 * @brief Implements SkyLine Structural Matrices. \ref structural "Structural Matrix"
 * ingroup structural
 */
class TPZSkylineStructMatrix : public TPZStructMatrix {
public:    
	
	TPZSkylineStructMatrix(TPZCompMesh *);
    
    TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);
	
	TPZSkylineStructMatrix(const TPZSkylineStructMatrix &cp);
    
    ~TPZSkylineStructMatrix();
	
    virtual TPZMatrix * Create();
	
    virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
    
    /** @brief Adapt the skyline for a range of equations */
    void FilterSkyline(TPZVec<int> &skyline);
	
public:
	
};

#endif //TPZSKYLINESTRUCTMATRIX_H
