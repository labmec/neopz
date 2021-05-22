/**
 * @file
 * @brief Contains TPZVerySparseMatrix class which implements a matrix whose nonzero elements are stored in binary tree.
 */
#ifndef TPZVERYSPARSEMATRIX_H
#define TPZVERYSPARSEMATRIX_H

#include "pzmatrix.h"
#include <iostream>
#include <map>
#include "Hash/TPZHash.h"

/** @ingroup matrix */

template <class TVar>
class TPZFMatrix;

template <class TVar>
class TPZFYsmpMatrix;

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
	
	TPZVerySparseMatrix(int64_t rows, int64_t cols) :
    TPZRegisterClassId(&TPZVerySparseMatrix::ClassId),
    TPZMatrix<TVar>(rows, cols)
	{
	}
	TPZVerySparseMatrix(int64_t rows, int64_t cols, TVar val) :
    TPZRegisterClassId(&TPZVerySparseMatrix::ClassId),
    TPZMatrix<TVar>(rows, cols)
	{
	}
	
	virtual ~TPZVerySparseMatrix();
	
	TPZVerySparseMatrix(const TPZVerySparseMatrix<TVar> &copy) :
    TPZRegisterClassId(&TPZVerySparseMatrix::ClassId),
    TPZMatrix<TVar>(copy), fExtraSparseData(copy.fExtraSparseData)
	{
	}
	
	TPZVerySparseMatrix(const TPZFMatrix<TVar> &cp);
	
	/** @brief Simetrizes copies the data of the matrix to make its data simetric */
	void Simetrize() override;
    
	/** @brief Put values checking bounds */
	int PutVal(const int64_t row, const int64_t col, const TVar &val) override;
	
	/** @brief Get values checking bounds */
	virtual const TVar GetVal(const int64_t row, const int64_t col) const override;
	
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
	 * @param row Row number.
	 * @param col Column number.
	 */
	virtual TVar &s(const int64_t row, const int64_t col) override
	{
#ifdef PZDEBUG
		if(row >= this->Rows() || row<0 || col >= this->Cols() || col<0)
		{
			this->Error("TPZFMatrix::operator() "," Index out of bounds");
			DebugStop();
		}
#endif
		return fExtraSparseData[std::pair<int64_t, int64_t>(row, col)];
	}
	inline TPZVerySparseMatrix<TVar>*NewMatrix() const override
  {return new TPZVerySparseMatrix<TVar>{};}
	CLONEDEF(TPZVerySparseMatrix)

	/** @brief Creates a copy from another TPZVerySparseMatrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override
  {                                                           
    auto *from = dynamic_cast<const TPZVerySparseMatrix<TVar> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 */
	virtual void MultAdd(const TPZFMatrix<TVar> & x, const TPZFMatrix<TVar> & y,
						 TPZFMatrix<TVar> & z, const TVar alpha = 1, const TVar beta = 0,
						 const int opt = 0) const override;
	
	/** @brief It makes *T the transpose of current matrix. */
	virtual void Transpose(TPZVerySparseMatrix<TVar>* T) const;
	virtual void Transpose(TPZMatrix<TVar>*const T) const  override {
        TPZMatrix<TVar>::Transpose(T);
    }
    
	/** @brief Saveable methods */
	public:
int ClassId() const override;

	void Write(TPZStream &buf, int withclassid) const override;
	void Read(TPZStream &buf, void *context) override;

	typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator MapBegin() const { return fExtraSparseData.begin(); }
	typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator MapEnd() const { return fExtraSparseData.end(); }
	
private:
	/** @brief Auxiliary functions only reading and writing a map as the third paremeter */
	void WriteMap(TPZStream &buf, int withclassid, const std::map<std::pair<int64_t, int64_t>, TVar> & TheMap) const;
	void ReadMap(TPZStream &buf, void *context, std::map<std::pair<int64_t, int64_t>, TVar> & TheMap);
	
protected:
	TVar *&Elem() override;
  const TVar *Elem() const override;
	inline int64_t Size() const override{
		return fExtraSparseData.size();
	}
    /** @brief Save elements different from zero, of Sparse matrix */
	std::map<std::pair<int64_t, int64_t>, TVar> fExtraSparseData;
    
};

template<class TVar>
int TPZVerySparseMatrix<TVar>::ClassId() const{
    return Hash("TPZVerySparseMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

#endif
