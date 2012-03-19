/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a symmetric sparse matrix.
 */
/******************************************************************************
 *
 * Class definition:    TPZSYsmpMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on symmetric sparse matrices stored
 *                      in the (old) Yale Sparse Matrix Package format.
 *
 * Solvers:             SOR
 *
 *****************************************************************************/

#ifndef SYSMPMATH
#define SYSMPMATH

#include "pzmatrix.h"
template<class TVar>
class TPZFMatrix;

/**
 * Purpose:  Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
 /**
  * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZSYsmpMatrix : public TPZMatrix<TVar>{
	
	public :
	
    TPZSYsmpMatrix( int rows, int cols );
	
    TPZSYsmpMatrix(const TPZSYsmpMatrix<TVar> &cp) : TPZMatrix<TVar>(cp)
    {
		int fjasize = fIA[this->Rows()];
		fIA = new int[this->Rows()+1];
		fDiag = new TVar[this->Rows()];
		fJA = new int[fjasize];
		fA = new TVar[fjasize];
		memcpy(fIA,cp.fIA,(this->Rows()+1)*sizeof(int));
		memcpy(fJA,cp.fJA,fjasize*sizeof(int));
		memcpy(fDiag,cp.fDiag,this->Rows()*sizeof(TVar));
		memcpy(fA,cp.fA,fjasize*sizeof(TVar));
		
    }
    
    CLONEDEF(TPZSYsmpMatrix)
	
	virtual ~TPZSYsmpMatrix();
	
	
	/// Get the matrix entry at (row,col) without bound checking
	virtual const TVar &GetVal(const int row,const int col ) const;
	
	/// computes z = beta * y + alpha * opt(this)*x \n
	///          z and x cannot overlap in memory
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 ) const ;
	
	/// Pass the data to the class.
	virtual void SetData( int *const IA, int *const JA, TVar *const A );
	
	/// Print the matrix along with a identification title
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const;
	
	void SolveSOR(int &numiterations,const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
				  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,
				  const TVar overrelax, TVar &tol,
				  const int FromCurrent = 0,const int direction = 1 ) ;
	
	
private:
	
	void ComputeDiagonal();
	
	int  *fIA;
	int  *fJA;
	TVar *fA;
	
	
	TVar *fDiag;
};

template<class TVar>
inline void TPZSYsmpMatrix<TVar>::SetData( int *IA, int *JA, TVar *A ) {
	// Pass the data to the class.
	fIA = IA;
	fJA = JA;
	fA  =  A;
	ComputeDiagonal();
}

#endif
