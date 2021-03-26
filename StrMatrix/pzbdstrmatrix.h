/**
 * @file
 * @brief Contains the TPZBlockDiagonalStructMatrix class which implements Block Diagonal Structural Matrices.
 */

#ifndef TPZBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZBLOCKDIAGONALSTRUCTMATRIX_H

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzstrmatrix.h"

#include "pzcmesh.h"
#include "pzvec.h"
#include "pzblockdiag.h"

/**
 * @brief Implements Block Diagonal Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZBlockDiagonalStructMatrix : public TPZStructMatrix {
public:
    
    int ClassId() const override;

    enum MBlockStructure {ENodeBased, EVertexBased, EElementBased};
    
    TPZBlockDiagonalStructMatrix(TPZCompMesh *);
    
    ~TPZBlockDiagonalStructMatrix();
    
    TPZBlockDiagonalStructMatrix(const TPZBlockDiagonalStructMatrix &copy) : 
    TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), TPZStructMatrix(copy),
    fBlockStructure(copy.fBlockStructure),fOverlap(copy.fOverlap)
    {
    }
    
    /** @brief Creates a sparse blockdiagonal matrix, overlapping should be assumed */
    virtual TPZMatrix<STATE> * Create() override;
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
    virtual TPZStructMatrix * Clone() override;
    
public:
    
    void AssembleBlockDiagonal(TPZBlockDiagonal<STATE> & block);
private:
    TPZBlockDiagonalStructMatrix();
    
    void BlockSizes(TPZVec < int > & blocksizes);
    
    friend TPZPersistenceManager;
    
    MBlockStructure fBlockStructure;
    int fOverlap;
    
    
};

#endif //TPZBLOCKDIAGONALSTRUCTMATRIX_H
