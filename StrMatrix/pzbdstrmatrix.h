/**
 * @file
 * @brief Contains the TPZBlockDiagonalStructMatrix class which implements Block Diagonal Structural Matrices.
 */

#ifndef TPZBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZBLOCKDIAGONALSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

template<class TVar>
class TPZBlockDiagonal;

/**
 * @brief Implements a block diagonal structural matrix.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZBlockDiagonalStructMatrix : public TPZStructMatrixT<TVar>,
                                     public TPar{
public:
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    
    enum MBlockStructure {ENodeBased, EVertexBased, EElementBased};   
    
    /** @brief Creates a sparse blockdiagonal matrix, overlapping should be assumed */
    TPZMatrix<TVar> * Create() override;
    
    TPZStructMatrix * Clone() override;

    void EndCreateAssemble(TPZBaseMatrix*) override;
    
    

    //@{
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
    
    void AssembleBlockDiagonal(TPZBlockDiagonal<TVar> & block);
private:
    
    void BlockSizes(TPZVec < int > & blocksizes);
    
    MBlockStructure fBlockStructure{EVertexBased};
    
    int fOverlap{0};
    
    friend TPZPersistenceManager;
};

#endif //TPZBLOCKDIAGONALSTRUCTMATRIX_H
