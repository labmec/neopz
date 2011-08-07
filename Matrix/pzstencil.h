/**
 * @file
 * @brief Contains TPZStencilMatrix class which implements a sparse matrix defined by a stencil.
 */
/******************************************************************************
 *
 * Class definition:    TPZStencilMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on sparse matrices stored by
 *                      stencils
 *
 * Solvers:             SOR
 *                      SSOR
 *
 *****************************************************************************/

#ifndef STENMATH
#define STENMATH

#include "pzmatrix.h"

class TPZFMatrix;

/**
 @brief Implements a sparse matrix defined by a stencil. \ref matrix "Matrix"
 @ingroup matrix
 */
 /**
  * Purpose:  Define operations on sparse matrices stored by stencils
  */
class TPZStencilMatrix : public TPZMatrix {
	
	public :
	
	/** @brief sets up the StencilMatrix based on the stencil */
    TPZStencilMatrix( int rows, int cols );
	
	virtual ~TPZStencilMatrix();
	
	/** @brief Returns the rows of this */
	virtual int Rows() const;
	/** @brief Returns the columns of this */
	virtual int Cols() const;
	
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const REAL &GetVal(const int row,const int col ) const;
	
	/** @brief computes \f$ z = beta * y + alpha * opt(this)*x \f$ */
	/**          z and x cannot overlap in memory */
	virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						 const REAL alpha=1., const REAL beta = 0., const int opt = 0 , const int stride = 1) const;
	
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
	
	void SolveSOR( int &numiterations,const TPZFMatrix &rhs, TPZFMatrix &x,
				  TPZFMatrix *residual, TPZFMatrix &scratch,
				  const REAL overrelax, REAL &tol,
				  const int FromCurrent = 0,const int direction = 1 );
	
	/** @brief initiates Stencil number "stencilnumber" with the data */
	void SetStencil( int stencilnumber, int inc, int *IA, REAL *A );
	
	/** @brief associates the given stencil number with each row */
	void SetNodeStencils( int *stencilnumber );
	
private:
	
	void IncreaseStencilPointers( int stencilnumber );
	
	void ComputeDiagonal();
	
	struct MPStencil {
		int fNumberOfItems;
		int fInc;
		int *fIA;
		REAL *fA;
		
		MPStencil( int inc, int *IA, REAL *A );
		~MPStencil();
		
	} **fMystencils;
	
	int   fNumberOfStencilPointers;
	int   fRows, fCols;
	int  *fStencilNumbers;
	
	REAL *fDiag;
	int   fSolver;
	int   fMaxIterations;
	
	REAL fSORRelaxation;
};

inline int TPZStencilMatrix::Rows() const{
	return fRows;
}

inline int TPZStencilMatrix::Cols() const{
	return fCols;
}

#endif
