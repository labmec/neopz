/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a symmetric sparse matrix. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */

#ifndef SYSMPMATH
#define SYSMPMATH

#include "pzmatrix.h"
template<class TVar>
class TPZFMatrix;

 /**
  * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZSYsmpMatrix : public TPZMatrix<TVar>{
	
	public :
    /** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrix();
	/** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrix(const long rows, const long cols );
	/** @brief Copy constructor */
    TPZSYsmpMatrix(const TPZSYsmpMatrix<TVar> &cp) : TPZMatrix<TVar>(cp), fIA(cp.fIA), fJA(cp.fJA), fA(cp.fA), fDiag(cp.fDiag)
    {
		
    }
    
    CLONEDEF(TPZSYsmpMatrix)
	/** @brief Destructor */
	virtual ~TPZSYsmpMatrix();
    
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(long nrow, long ncol, int symmetric);
	
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const TVar &GetVal(const long row, const long col ) const;
	
	/** @brief Computes z = beta * y + alpha * opt(this)*x */
	/** @note z and x cannot overlap in memory */
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const ;
	
	/** @brief Sets data to the class */
	virtual void SetData(const TPZVec<long> &IA,const TPZVec<long> &JA, const TPZVec<TVar> &A );
	
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const;
	
	
private:
	
	void ComputeDiagonal();
	
	TPZVec<long>  fIA;
	TPZVec<long>  fJA;
	TPZVec<TVar> fA;
	
	
	TPZVec<TVar> fDiag;
};

template<class TVar>
inline void TPZSYsmpMatrix<TVar>::SetData(const TPZVec<long> &IA,const TPZVec<long> &JA, const TPZVec<TVar> &A )
{
	// Pass the data to the class.
	fIA = IA;
	fJA = JA;
	fA  =  A;
	ComputeDiagonal();
}

#endif
