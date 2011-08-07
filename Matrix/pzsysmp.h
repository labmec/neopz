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

class TPZFMatrix;

/**
 * Purpose:  Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
 /**
  * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
  * @ingroup matrix
  */
class TPZSYsmpMatrix : public TPZMatrix {
	
	public :
	
    TPZSYsmpMatrix( int rows, int cols );
	
    TPZSYsmpMatrix(const TPZSYsmpMatrix &cp) : TPZMatrix(cp)
    {
		int fjasize = fIA[Rows()];
		fIA = new int[Rows()+1];
		fDiag = new REAL[Rows()];
		fJA = new int[fjasize];
		fA = new REAL[fjasize];
		memcpy(fIA,cp.fIA,(Rows()+1)*sizeof(int));
		memcpy(fJA,cp.fJA,fjasize*sizeof(int));
		memcpy(fDiag,cp.fDiag,Rows()*sizeof(REAL));
		memcpy(fA,cp.fA,fjasize*sizeof(REAL));
		
    }
    
    CLONEDEF(TPZSYsmpMatrix)
	
	virtual ~TPZSYsmpMatrix();
	
	
	/// Get the matrix entry at (row,col) without bound checking
	virtual const REAL &GetVal(const int row,const int col ) const;
	
	/// computes z = beta * y + alpha * opt(this)*x \n
	///          z and x cannot overlap in memory
	virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						 const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 ) const ;
	
	/// Pass the data to the class.
	virtual void SetData( int *const IA, int *const JA, REAL *const A );
	
	/// Print the matrix along with a identification title
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const;
	
	void SolveSOR(int &numiterations,const TPZFMatrix &rhs, TPZFMatrix &x,
				  TPZFMatrix *residual, TPZFMatrix &scratch,
				  const REAL overrelax, REAL &tol,
				  const int FromCurrent = 0,const int direction = 1 ) ;
	
	
private:
	
	void ComputeDiagonal();
	
	int  *fIA;
	int  *fJA;
	REAL *fA;
	
	
	REAL *fDiag;
};


inline void TPZSYsmpMatrix::SetData( int *IA, int *JA, REAL *A ) {
	// Pass the data to the class.
	fIA = IA;
	fJA = JA;
	fA  =  A;
	ComputeDiagonal();
}

#endif
