/**
 * @file
 * @brief Contains TPZVerySparseMatrix class which implements a matrix whose nonzero elements are stored in binary tree.
 */
#ifndef TPZVERYSPARSEMATRIX_H
#define TPZVERYSPARSEMATRIX_H

#include <iostream>
#include <map>

#include "pzmatrix.h"
#include "pzfmatrix.h"
template<class TVar>
class TPZFYsmpMatrix;

/** @ingroup matrix */
#define TPZVERYSPARSEMATRIX_ID 28291001;

/**
 * @author Agnaldo Monteiro Farias <agnaldo@labmec.fec.unicamp.br>
 * @brief Implements a matrix whose nonzero elements are stored in binary tree. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZVerySparseMatrix: public TPZMatrix<TVar>
{
public:
	
	friend class TPZFYsmpMatrix<TVar>;
	
	TPZVerySparseMatrix();
	
	TPZVerySparseMatrix(long rows, long cols) :
    TPZMatrix<TVar>(rows, cols)
	{
	}
	TPZVerySparseMatrix(long rows, long cols, TVar val) :
    TPZMatrix<TVar>(rows, cols)
	{
	}
	
	virtual ~TPZVerySparseMatrix();
	
	TPZVerySparseMatrix(const TPZVerySparseMatrix<TVar> &copy) :
    TPZMatrix<TVar>(copy), fExtraSparseData(copy.fExtraSparseData)
	{
	}
	
	TPZVerySparseMatrix(const TPZFMatrix<TVar> &cp);
	
	/** @brief Simetrizes copies the data of the matrix to make its data simetric */
	void Simetrize();
    
	/** @brief Put values checking bounds */
	int PutVal(const long row, const long col, const TVar &val);
	
	/** @brief Get values checking bounds */
	virtual const TVar &GetVal(const long row, const long col) const;
	
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
	 * @param row Row number.
	 * @param col Column number.
	 */
	virtual TVar &s(const long row, const long col)
	{
#ifdef DEBUG
		if(row >= this->Rows() || row<0 || col >= this->Cols() || col<0)
		{
			this->Error("TPZFMatrix::operator() "," Index out of bounds");
			DebugStop();
			return this->gZero;
		}
#endif
		return fExtraSparseData[std::pair<long, long>(row, col)];
	}
	
	CLONEDEF(TPZVerySparseMatrix)
	
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
	 */
	virtual void MultAdd(const TPZFMatrix<TVar> & x, const TPZFMatrix<TVar> & y,
						 TPZFMatrix<TVar> & z, const TVar alpha = 1, const TVar beta = 0,
						 const int opt = 0, const int stride = 1) const;
	
	/** @brief It makes *T the transpose of current matrix. */
	virtual void Transpose(TPZVerySparseMatrix<TVar>* T) const;
	virtual void Transpose(TPZMatrix<TVar>*const T) const {
        TPZMatrix<TVar>::Transpose(T);
    }
    
	/** @brief Saveable methods */
	int ClassId() const
	{
		return TPZVERYSPARSEMATRIX_ID;
	}
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);

	typename std::map <std::pair<long, long>, TVar>::const_iterator MapBegin() const { return fExtraSparseData.begin(); }
	typename std::map <std::pair<long, long>, TVar>::const_iterator MapEnd() const { return fExtraSparseData.end(); }
	
private:
	/** @brief Auxiliary functions only reading and writing a map as the third paremeter */
	void WriteMap(TPZStream &buf, int withclassid, std::map<std::pair<long, long>, TVar> & TheMap);
	void ReadMap(TPZStream &buf, void *context, std::map<std::pair<long, long>, TVar> & TheMap);
	
protected:
    
    /** @brief Save elements different from zero, of Sparse matrix */
	std::map<std::pair<long, long>, TVar> fExtraSparseData;
    
};

#endif
