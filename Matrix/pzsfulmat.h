/**
 * @file
 * @brief Contains TPZSFMatrix class which implements a symmetric full matrix.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tsfulmat.hh
//
// Class:  TPZSFMatrix
//
// Obs.:   Implementa matrizes cheias simetricas.
//
// Versao: 04 / 1996.
//


#ifndef _TSFULMATHH_
#define _TSFULMATHH_

#include "pzfmatrix.h"
#include "pzmatrix.h"
//#include "pztempmat.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

class TPZFMatrix;

/**
 @brief Implements a symmetric full matrix. \ref matrix "Matrix"
 @ingroup matrix
 */
class TPZSFMatrix : public TPZMatrix {
	//  friend class TSimMatrix;
	
public:
	TPZSFMatrix () : TPZMatrix( 0,0 )  { fElem = NULL; }
	TPZSFMatrix (const int dim );
	TPZSFMatrix (const TPZSFMatrix & );
	// Usa o maior bloco quadrado possivel, comecado em (0,0).
	// E inicializa com a parte triangular inferior do bloco.
	TPZSFMatrix (const TPZMatrix & );
	
	CLONEDEF(TPZSFMatrix)
	
	~TPZSFMatrix();
	
	int PutVal(const int row,const int col,const REAL &value );
	const REAL &GetVal(const int row,const int col ) const;
	
	// Peforms the product (*this)T x D x (*this).
	//TPZSFMatrix  InnerProd(TPZSFMatrix & D);
	
	/// Operators with Full simmetric matrices.
	// @{
	TPZSFMatrix &operator= (const TPZSFMatrix &A );
	TPZSFMatrix operator+  (const TPZSFMatrix &A ) const;
	TPZSFMatrix operator-  (const TPZSFMatrix &A ) const;
	TPZSFMatrix &operator+=(const TPZSFMatrix &A );
	TPZSFMatrix &operator-=(const TPZSFMatrix &A );
	// @}
	
	/// Operators with generic matrices.
	// @{
	TPZSFMatrix &operator= (const TPZMatrix &A );
	TPZSFMatrix operator+  (const TPZMatrix &A ) const;
	TPZSFMatrix operator-  (const TPZMatrix &A ) const;
	TPZSFMatrix &operator+=(const TPZMatrix &A );
	TPZSFMatrix &operator-=(const TPZMatrix &A );
	// @}
	
	//TTempMat<TPZFMatrix> operator+(const TPZMatrix &A ) const {return TPZMatrix::operator+(A);}
	//TTempMat<TPZFMatrix> operator-(const TPZMatrix &A ) const {return TPZMatrix::operator-(A);}
	//TTempMat<TPZFMatrix> operator*( TPZFMatrix &A ) const {return TPZMatrix::operator*(A);}
	
	/// Operators with numeric values.
	// @{
	TPZSFMatrix &operator= (const REAL val );
	TPZSFMatrix operator+  (const REAL val ) const;
	TPZSFMatrix operator-  (const REAL val ) const { return operator+( -val ); }
	TPZSFMatrix operator*  (const REAL val ) const;
	TPZSFMatrix &operator+=(const REAL val );
	TPZSFMatrix &operator-=(const REAL val )  { return operator+=( -val ); }
	TPZSFMatrix &operator*=(const REAL val );
	// @}
	
	TPZSFMatrix operator-() const  { return operator*( -1.0 ); }
	
	/// Resize the array but keeps its entirety.
	int Resize(const int newDim, const int );
	
	// Resize the array and resets ist entirety.
	int Redim(const int newRows ,const int);
	
	int Redim(const int newDim) {return Redim(newDim,newDim);}
	
	// Resets all elements.
	int Zero();
	
	
	/*** To solve linear systems ***/
	// @{
	virtual int Decompose_Cholesky();
	virtual int Decompose_LDLt();
	virtual int Decompose_Cholesky(std::list<int> &singular);
	virtual int Decompose_LDLt(std::list<int> &singular);
	
	virtual int Subst_Forward  ( TPZFMatrix *B ) const;
	virtual int Subst_Backward ( TPZFMatrix *B ) const;
	virtual int Subst_LForward ( TPZFMatrix *B ) const;
	virtual int Subst_LBackward( TPZFMatrix *B ) const;
	virtual int Subst_Diag     ( TPZFMatrix *B ) const;
	// @}
	
#ifdef OOPARLIB
	
	virtual long GetClassID() const  { return TSFMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual char *ClassName() const   { return( "TPZSFMatrix"); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const;
	
#endif
	
private:
	
	int Size() const { return (Dim() * (Dim()+1)) >> 1; }
	
	//static  int Error(const char *msg1,const char *msg2="" );
	int Clear();
	
	REAL   *fElem;
};

/**************/
/*** PutVal ***/
inline int
TPZSFMatrix::PutVal(const int row,const int col,const REAL &value )
{
	int locrow = row, loccol = col;
	if ( locrow < loccol )
		Swap( &locrow, &loccol );
	
	fDecomposed = 0;
	fElem[ ((locrow*(locrow+1)) >> 1) + loccol ] = value;
	return( 1 );
}

/**************/
/*** GetVal ***/
inline const REAL &
TPZSFMatrix::GetVal(const int row,const int col ) const
{
	int locrow(row),loccol(col);
	if ( locrow < loccol )
		Swap( &locrow, &loccol );
	
	return( fElem[ ((locrow*(locrow+1)) >> 1) + loccol ] );
}

/**************************/
/*** Operadores Globais ***/

/** @brief Increments value for all entries of the A matrix */
inline TPZSFMatrix
operator+( const REAL value, const TPZSFMatrix &A )
{
	return( A + value );
}

/** @brief Decrements value for all entries of the A matrix */
inline TPZSFMatrix
operator-( const REAL value,const  TPZSFMatrix &A )
{
	return( A - value );
}

/** @brief Implements the scalar product value*A */
inline TPZSFMatrix
operator*( const REAL value, const TPZSFMatrix &A )
{
	return( A * value );
}

#endif

