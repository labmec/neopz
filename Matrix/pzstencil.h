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

template<class TVar>
class TPZFMatrix;

/**
 @brief Implements a sparse matrix defined by a stencil. \ref matrix "Matrix"
 @ingroup matrix
 */
 /**
  * Purpose:  Define operations on sparse matrices stored by stencils
  */
template<class TVar>
class TPZStencilMatrix : public TPZMatrix<TVar> {
	
	public :
	
	/** @brief sets up the StencilMatrix based on the stencil */
    TPZStencilMatrix( int rows, int cols );
	
	virtual ~TPZStencilMatrix();
	
	/** @brief Returns the rows of this */
	virtual int Rows() const;
	/** @brief Returns the columns of this */
	virtual int Cols() const;
	
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const TVar &GetVal(const int row,const int col ) const;
	
	/** @brief computes \f$ z = beta * y + alpha * opt(this)*x \f$ */
	/**          z and x cannot overlap in memory */
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1., const TVar beta = 0., const int opt = 0 , const int stride = 1) const;
	
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
	
	void SolveSOR( int &numiterations,const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
				  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,
				  const TVar overrelax, TVar &tol,
				  const int FromCurrent = 0,const int direction = 1 );
	
	/** @brief initiates Stencil number "stencilnumber" with the data */
	void SetStencil( int stencilnumber, int inc, int *IA, TVar *A );
	
	/** @brief associates the given stencil number with each row */
	void SetNodeStencils( int *stencilnumber );
	
private:
	
	void IncreaseStencilPointers( int stencilnumber );
	
	void ComputeDiagonal();
	
	struct MPStencil {
		int fNumberOfItems;
		int fInc;
		int *fIA;
		TVar *fA;
		
		MPStencil( int inc, int *IA, TVar *A );
		~MPStencil();
		
	} **fMystencils;
	
	int   fNumberOfStencilPointers;
	int   fRows, fCols;
	int  *fStencilNumbers;
	
	TVar *fDiag;
	int   fSolver;
	int   fMaxIterations;
	
	TVar fSORRelaxation;
};

template<class TVar>
inline int TPZStencilMatrix<TVar>::Rows() const{
	return fRows;
}

template<class TVar>
inline int TPZStencilMatrix<TVar>::Cols() const{
	return fCols;
}

#endif
