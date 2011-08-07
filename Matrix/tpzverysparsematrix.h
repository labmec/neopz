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
class TPZFYsmpMatrix;

/** @ingroup matrix */
#define TPZVERYSPARSEMATRIX_ID 28291001;

/**
 @author Agnaldo Monteiro Farias <agnaldo@labmec.fec.unicamp.br>
 @brief Implements a matrix whose nonzero elements are stored in binary tree. \ref matrix "Matrix"
 @ingroup matrix
 */
class TPZVerySparseMatrix: public TPZMatrix
{
public:
	
	friend class TPZFYsmpMatrix;
	
	TPZVerySparseMatrix();
	
	TPZVerySparseMatrix(int rows, int cols) :
    TPZMatrix(rows, cols)
	{
	}
	TPZVerySparseMatrix(int rows, int cols, REAL val) :
    TPZMatrix(rows, cols)
	{
	}
	
	virtual ~TPZVerySparseMatrix();
	
	TPZVerySparseMatrix(const TPZVerySparseMatrix &copy) :
    TPZMatrix(copy), fExtraSparseData(copy.fExtraSparseData)
	{
	}
	
	TPZVerySparseMatrix(const TPZFMatrix &cp);
	
	/**
	 * @brief Put values checking bounds\n
	 */
	virtual int PutVal(const int row, const int col, const REAL &val);
	
	/**
	 * @brief Get values checking bounds\n
	 */
	virtual const REAL &GetVal(const int row, const int col) const;
	
	/**
	 @brief The operators check on the bounds if the DEBUG variable is defined
	 @param row Row number.
	 @param col Column number.
	 */
	virtual REAL &s(const int row, const int col)
	{
#ifdef DEBUG
		if(row >= Rows() || row<0 || col >= Cols() || col<0)
		{
			Error("TPZFMatrix::operator() "," Index out of bounds");
			DebugStop();
			gZero = 0.;
			return gZero;
		}
#endif
		return fExtraSparseData[std::pair<int, int>(row, col)];
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
	virtual void MultAdd(const TPZFMatrix & x, const TPZFMatrix & y,
						 TPZFMatrix & z, const REAL alpha = 1., const REAL beta = 0.,
						 const int opt = 0, const int stride = 1) const;
	
	/**
	 * @brief Saveable methods
	 */
	int ClassId() const
	{
		return TPZVERYSPARSEMATRIX_ID;
	}
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
private:
	/**
	 * @brief Auxiliary functions only reading and writing a map as the third paremeter
	 */
	void WriteMap(TPZStream &buf, int withclassid, std::map<std::pair<int, int>, REAL> & TheMap);
	void ReadMap(TPZStream &buf, void *context, std::map<std::pair<int, int>, REAL> & TheMap);
	
	
	
protected:
	
	/**
	 @brief Save elements different from zero, of Sparse matrix
	 */
	std::map<std::pair<int, int>, REAL> fExtraSparseData;
	
};

#endif
