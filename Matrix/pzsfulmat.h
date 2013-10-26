/**
 * @file
 * @brief Contains TPZSFMatrix class which implements a symmetric full matrix.
 */

#ifndef _TSFULMATHH_
#define _TSFULMATHH_

#include "pzfmatrix.h"
#include "pzmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif
template<class TVar>
class TPZFMatrix;

/**
 * @brief Implements a symmetric full matrix. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSFMatrix : public TPZMatrix<TVar> {
	
public:
	TPZSFMatrix () : TPZMatrix<TVar>( 0,0 )  { fElem = NULL; }
	TPZSFMatrix (const long dim );
	TPZSFMatrix (const TPZSFMatrix<TVar> & );
	// Usa o maior bloco quadrado possivel, comecado em (0,0).
	// E inicializa com a parte triangular inferior do bloco.
	TPZSFMatrix (const TPZMatrix<TVar> & );
	
	CLONEDEF(TPZSFMatrix)
	
	~TPZSFMatrix();
	
	int PutVal(const long row,const long col,const TVar &value );
	const TVar &GetVal(const long row,const long col ) const;
	
	/**
	 * @name Operators with Full simmetric matrices.
	 * @{
	 */
	TPZSFMatrix &operator= (const TPZSFMatrix<TVar> &A );
	TPZSFMatrix operator+  (const TPZSFMatrix<TVar> &A ) const;
	TPZSFMatrix operator-  (const TPZSFMatrix<TVar> &A ) const;
	TPZSFMatrix &operator+=(const TPZSFMatrix<TVar> &A );
	TPZSFMatrix &operator-=(const TPZSFMatrix<TVar> &A );
	/** @} */
	
	/**
	 * @name Operators with generic matrices.
	 * @{
	 */
	TPZSFMatrix &operator= (const TPZMatrix<TVar> &A );
	TPZSFMatrix operator+  (const TPZMatrix<TVar> &A ) const;
	TPZSFMatrix operator-  (const TPZMatrix<TVar> &A ) const;
	TPZSFMatrix &operator+=(const TPZMatrix<TVar> &A );
	TPZSFMatrix &operator-=(const TPZMatrix<TVar> &A );
	/** @} */

	/**
	 * @name Operators with numeric values.
	 * @{
	 */
	TPZSFMatrix &operator= (const TVar val );
	TPZSFMatrix operator+  (const TVar val ) const;
	TPZSFMatrix operator-  (const TVar val ) const { return operator+( -val ); }
	TPZSFMatrix operator*  (const TVar val ) const;
	TPZSFMatrix &operator+=(const TVar val );
	TPZSFMatrix &operator-=(const TVar val )  { return operator+=( -val ); }
	TPZSFMatrix &operator*=(const TVar val );
	/** @} */
	
	TPZSFMatrix operator-() const  { return operator*( -1.0 ); }
	
	/** @brief Resize the array but keeps its entirety. */
	int Resize(const long newDim, const long );
	
	/** @brief Resize the array and resets ist entirety. */
	int Redim(const long newRows ,const long);
	
	int Redim(const long newDim) {return Redim(newDim,newDim);}
	
	/** @brief Resets all elements. */
	int Zero();
	
	
	/***
	 * @name To solve linear systems
	 * @{
	 */
	virtual int Decompose_Cholesky();
	virtual int Decompose_LDLt();
	virtual int Decompose_Cholesky(std::list<long> &singular);
	virtual int Decompose_LDLt(std::list<long> &singular);
	
	virtual int Subst_Forward  ( TPZFMatrix<TVar>  *B ) const;
	virtual int Subst_Backward ( TPZFMatrix<TVar>  *B ) const;
	virtual int Subst_LForward ( TPZFMatrix<TVar>  *B ) const;
	virtual int Subst_LBackward( TPZFMatrix<TVar>  *B ) const;
	virtual int Subst_Diag     ( TPZFMatrix<TVar>  *B ) const;
	/** @} */
	
#ifdef OOPARLIB
	
	virtual long GetClassID() const  { return TSFMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZSFMatrix"); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const;
	
#endif
	
private:
	
	long Size() const { return (this->Dim() * (this->Dim()+1)) >> 1; }
	
	int Clear();
	
	TVar   *fElem;
};

/**************/
/*** PutVal ***/
template<class TVar>
inline int
TPZSFMatrix<TVar> ::PutVal(const long row,const long col,const TVar &value )
{
	long locrow = row, loccol = col;
	if ( locrow < loccol )
		this->Swap( &locrow, &loccol );
	
	this->fDecomposed = 0;
	this->fElem[ ((locrow*(locrow+1)) >> 1) + loccol ] = value;
	return( 1 );
}

/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &
TPZSFMatrix<TVar> ::GetVal(const long row,const long col ) const
{
	long locrow(row),loccol(col);
	if ( locrow < loccol )
		this->Swap( &locrow, &loccol );
	
	return( fElem[ ((locrow*(locrow+1)) >> 1) + loccol ] );
}

/**************************/
/*** @name Operadores Globais ***/
/** @{ */

/** @brief Increments value for all entries of the A matrix */
template<class TVar>
inline TPZSFMatrix<TVar> 
operator+( const TVar value, const TPZSFMatrix<TVar>  &A )
{
	return( A + value );
}

/** @brief Decrements value for all entries of the A matrix */
template<class TVar>
inline TPZSFMatrix<TVar> 
operator-( const TVar value,const  TPZSFMatrix<TVar>  &A )
{
	return( A - value );
}

/** @brief Implements the scalar product value*A */
template<class TVar>
inline TPZSFMatrix<TVar> 
operator*( const TVar value, const TPZSFMatrix<TVar>  &A )
{
	return( A * value );
}

/** @} */

#endif
