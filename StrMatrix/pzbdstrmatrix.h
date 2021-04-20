/**
 * @file
 * @brief Contains the TPZBlockDiagonalStructMatrix class which implements Block Diagonal Structural Matrices.
 */

#ifndef TPZBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZBLOCKDIAGONALSTRUCTMATRIX_H

#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"
template<class TVar>
class TPZBlockDiagonal;

/**
 * @brief Implements Block Diagonal Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZBlockDiagonalStructMatrix :
    public TPZStructMatrix,
    public TPZStructMatrixOR<STATE>{
public:

    enum MBlockStructure {ENodeBased, EVertexBased, EElementBased};
    
    TPZBlockDiagonalStructMatrix(TPZCompMesh *);
    TPZBlockDiagonalStructMatrix(TPZAutoPointer<TPZCompMesh>);
    
    ~TPZBlockDiagonalStructMatrix();
    
    TPZBlockDiagonalStructMatrix(const TPZBlockDiagonalStructMatrix &copy) : 
    TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), TPZStructMatrix(copy),
    fBlockStructure(copy.fBlockStructure),fOverlap(copy.fOverlap)
    {
    }
    
    /** @brief Creates a sparse blockdiagonal matrix, overlapping should be assumed */
    TPZMatrix<STATE> * Create() override;
    
    TPZStructMatrix * Clone() override;
    
    TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
    

    //@{
    //!Read and Write methods
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
    
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
