/**
 * @file
 * @brief Contains TPZSparseBlockDiagonal class which implements a block diagonal matrix where the blocks are not contiguous.
 */
//
// C++ Interface: tpzsparseblockdiagonal
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZSPARSEBLOCKDIAGONAL_H
#define TPZSPARSEBLOCKDIAGONAL_H

#include "pzblockdiag.h"

/**
 @brief Implements a block diagonal matrix where the blocks are not contiguous. \ref matrix "Matrix"
 @ingroup matrix
 @author Philippe R. B. Devloo
 @since 2004
 */
template<class TVar>
class TPZSparseBlockDiagonal : public TPZBlockDiagonal<TVar>
{
public:
    TPZSparseBlockDiagonal();
    TPZSparseBlockDiagonal(TPZVec<int> &blockgraph, TPZVec<int> &blockgraphindex,int rows);
    
    TPZSparseBlockDiagonal(TPZVec<int> &blockgraph, TPZVec<int> &blockgraphindex,int rows, int color, TPZVec<int> &colors);
	
    ~TPZSparseBlockDiagonal();
	
    const TVar& Get(const int row, const int col) const;
    const TVar& GetVal(const int row, const int col) const;
    int Put(const int row, const int col, const TVar& value);
    int PutVal(const int row, const int col, const TVar& value);
    TVar& operator ( )(const int row, const int col);
    virtual int Substitution(TPZFMatrix<TVar>* B) const;
    virtual TVar& s(const int row, const int col);
    virtual void Print(const char* message, std::ostream& out, MatrixOutputFormat=EFormatted) const;
    void AddBlock(int i, TPZFMatrix<TVar>& block);
    void BuildFromMatrix(TPZMatrix<TVar>& matrix);
    void GetBlock(int i, TPZFMatrix<TVar>& block);
    void MultAdd(const TPZFMatrix<TVar>& x, const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z, const TVar alpha, const TVar beta, const int opt, const int stride) const;
    void FindBlockIndex(int glob, int &block, int &blockind) const;
	
	/**
	 * @brief Updates the values of the matrix based on the values of the matrix
	 */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
	
    
protected:
	/**
	 @brief Equation numbers for each block
	 */
    TPZVec<int> fBlock;
	/**
	 @brief Index to first element of each block
	 */
    TPZVec<int> fBlockIndex;
	
    void Scatter(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out, int stride) const;
    void Gather(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out, int stride) const;
};

#endif
