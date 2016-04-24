/**
 * @file
 * @brief Contains TPZSBMatrixLapack class which implements symmetric band matrices using resources from LAPACK.
 * 
 * @details While the class TPZSBMatrix also implements a symmetric band matrix, it does so in a manner not adequate
 * for straightforward use of LAPACK methods, since it stores the relevant diagonals in a column-major fashion, beginning
 * with the main diagonal. LAPACK requires the matrix to be stored in a row-major fashion, beginning with the
 * diagonal in which j - i = fBand
 */

#ifndef TSBNDMATH
#define TSBNDMATH

#include "pzmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

template<class TVar>
class TPZFMatrix;

/**
 * @brief Implements symmetric band matrices. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSBMatrixLapack : public TPZMatrix<TVar>
{
public:
  TPZSBMatrixLapack();
	TPZSBMatrixLapack(const long dim,const long band );
  TPZSBMatrixLapack(const TPZSBMatrixLapack<TVar> &A );
	
	CLONEDEF(TPZSBMatrixLapack)
	
	~TPZSBMatrixLapack() { Clear(); }
	
	int    PutVal(const long row,const long col,const TVar& element );
	const TVar &GetVal(const long row,const long col ) const;
	
	/** @brief Computes z = beta * y + alpha * opt(this)*x */
	/** z and x cannot overlap in memory */
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 ) const;
	
	void Print(const char *name = NULL, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
	//friend std::ostream & operator<< <>(std::ostream& out,const TPZSBMatrixLapack<TVar>  &A); Leonardo removendo o '<>' antes do (std...
	template<class TT>friend std::ostream & operator<< (std::ostream& out,const TPZSBMatrixLapack<TT>  &A); 
	
	/// Operadores com matrizes SKY LINE.
	// @{
	TPZSBMatrixLapack &operator= (const TPZSBMatrixLapack<TVar> &A );
	TPZSBMatrixLapack operator+  (const TPZSBMatrixLapack<TVar> &A ) const;
	TPZSBMatrixLapack operator-  (const TPZSBMatrixLapack<TVar> &A ) const;
	TPZSBMatrixLapack &operator+=(const TPZSBMatrixLapack<TVar> &A );
	TPZSBMatrixLapack &operator-=(const TPZSBMatrixLapack<TVar> &A );
	// @}
	TPZSBMatrixLapack<TVar> operator*  (const TVar v ) const;
	TPZSBMatrixLapack<TVar> &operator+=(const TVar v );
	TPZSBMatrixLapack<TVar> &operator*=(const TVar v );
	
	TPZSBMatrixLapack<TVar> operator-() const { return operator*(-1.0); }
	
	/// Redimension the matrix keeping original elements.
	int Resize(const long newDim ,const long);
	
	/// Redimension the matrix and zeroes its elements
	int Redim(const long newDim) {return Redim(newDim,newDim);}
	int Redim(const long newRows ,const long newCols);
	
	/// Zeroes the elements of the matrix
	int Zero();
	
	long GetBand() const { return fBand; }
	int SetBand(const long newBand );
	
	
	/// To solve linear systems
	// @{
	int Decompose_Cholesky();  // Faz A = GGt.
	int Decompose_LDLt    ();  // Faz A = LDLt.
	int Decompose_Cholesky(std::list<long> &singular);  // Faz A = GGt.
	int Decompose_LDLt    (std::list<long> &singular);  // Faz A = LDLt.
	
	int Subst_Forward  ( TPZFMatrix<TVar> *b ) const;
	int Subst_Backward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LBackward( TPZFMatrix<TVar> *b ) const;
	int Subst_Diag     ( TPZFMatrix<TVar> *b ) const;
	// @}
	
	
private:
	
	long  Size() const    { return( this->Dim() * (fBand + 1) ); }
	int  PutZero();
	int  Clear();
	void Copy (const TPZSBMatrixLapack<TVar> & );
	
	
	TVar *fDiag;
	long  fBand;
};

#endif
