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
    TPZSparseBlockDiagonal(TPZVec<long> &blockgraph, TPZVec<long> &blockgraphindex,long rows);
    
    TPZSparseBlockDiagonal(TPZVec<long> &blockgraph, TPZVec<long> &blockgraphindex,long rows, int color, TPZVec<int> &colors);
	
    ~TPZSparseBlockDiagonal();
	
    const TVar& Get(const long row, const long col) const;
    const TVar& GetVal(const long row, const long col) const;
    int Put(const long row, const long col, const TVar& value);
    int PutVal(const long row, const long col, const TVar& value);
    TVar& operator ( )(const long row, const long col);
    virtual int Substitution(TPZFMatrix<TVar>* B) const;
    
    virtual TVar &s(const long row, const long col);
    
    virtual void Print(const char* message, std::ostream& out, const MatrixOutputFormat=EFormatted) const;
    void AddBlock(long i, TPZFMatrix<TVar>& block);
    void BuildFromMatrix(TPZMatrix<TVar>& matrix);
    void GetBlock(long i, TPZFMatrix<TVar>& block);
    void MultAdd(const TPZFMatrix<TVar>& x, const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z, const TVar alpha, const TVar beta, const int opt, const int stride) const;
    void FindBlockIndex(long glob, long &block, long &blockind) const;
	
	/** @brief Updates the values of the matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
	
    
protected:
	/** @brief Equation numbers for each block */
    TPZVec<long> fBlock;
	/** @brief Index to first element of each block */
    TPZVec<long> fBlockIndex;
	
    void Scatter(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out, int stride) const;
    void Gather(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out, int stride) const;
};

#endif
