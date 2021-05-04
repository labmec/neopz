/**
 * @file
 * @brief Contains TPZSparseBlockDiagonal class which implements a block diagonal matrix where the blocks are not contiguous.
 */
#ifndef TPZSPARSEBLOCKDIAGONAL_H
#define TPZSPARSEBLOCKDIAGONAL_H

#include "pzblockdiag.h"

/**
 * @brief Implements a block diagonal matrix where the blocks are not contiguous. \ref matrix "Matrix"
 * @ingroup matrix
 * @author Philippe R. B. Devloo
 * @since 2004
 */
template<class TVar>
class TPZSparseBlockDiagonal : public TPZBlockDiagonal<TVar>
{
public:
    TPZSparseBlockDiagonal();
    TPZSparseBlockDiagonal(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex,int64_t rows);
    
    TPZSparseBlockDiagonal(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex,int64_t rows, int color, TPZVec<int> &colors);
	
    ~TPZSparseBlockDiagonal();
	
    const TVar Get(const int64_t row, const int64_t col) const override;
    const TVar GetVal(const int64_t row, const int64_t col) const override;
    int Put(const int64_t row, const int64_t col, const TVar& value) override;
    int PutVal(const int64_t row, const int64_t col, const TVar& value) override;
    TVar& operator ( )(const int64_t row, const int64_t col);
    virtual int Substitution(TPZFMatrix<TVar>* B) const override;
    
    virtual TVar &s(const int64_t row, const int64_t col) override;
    
    virtual void Print(const char* message, std::ostream& out, const MatrixOutputFormat=EFormatted) const override;
    void AddBlock(int64_t i, TPZFMatrix<TVar>& block);
    void BuildFromMatrix(TPZMatrix<TVar>& matrix);
    void GetBlock(int64_t i, TPZFMatrix<TVar>& block);
    void MultAdd(const TPZFMatrix<TVar>& x, const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z, const TVar alpha, const TVar beta, const int opt) const override;
    void FindBlockIndex(int64_t glob, int64_t &block, int64_t &blockind) const;
	
	/** @brief Updates the values of the matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat) override;
	
    public:
int ClassId() const override;

protected:
	/** @brief Equation numbers for each block */
    TPZVec<int64_t> fBlock;
	/** @brief Index to first element of each block */
    TPZVec<int64_t> fBlockIndex;
	
    void Scatter(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out) const;
    void Gather(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out) const;
};

#endif
