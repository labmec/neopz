/**
 * @file
 * @brief Contains TPZFBMatrix class which defines a non symmetric banded matrix.
 */

#ifndef _TBNDMATHH_
#define _TBNDMATHH_

#include "pzmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

/**
 * @brief Defines a non symmetric banded matrix. \ref matrix "Matrix"
 * @ingroup matrix
 * @author MISAEL LUIS SANTANA MANDUJANO.
 * @since 12/1994
 */
template<class TVar>
class TPZFBMatrix : public TPZMatrix<TVar>
{
	
public:
	virtual int Substitution(TPZFMatrix<TVar> *B) const;
	/** @brief Simple constructor */
	TPZFBMatrix ();
	/**
     * @brief Constructor
     * @param dim Initial dimension of banded matrix
     * @param band_width Initial band width of banded matrix
	 */
	TPZFBMatrix (const long dim,const long band_width = 0 );
	/** @brief Copy constructor */
	TPZFBMatrix (const TPZFBMatrix<TVar> & );
	
	CLONEDEF(TPZFBMatrix)
	/** @brief Simple destructor */
	~TPZFBMatrix();

	int    Put(const long row,const long col,const TVar& value );
	const TVar &Get(const long row,const long col ) const;
	
	TVar &operator()(const long row, const long col);
	virtual TVar &s(const long row, const long col);

	inline int    PutVal(const long row,const long col,const TVar& value );
	inline const TVar &GetVal(const long row,const long col ) const;
	
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha=1,const TVar beta = 0,const int opt = 0,const int stride = 1 ) const;
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	
	// Peforms the product (*this)T x D x (*this).
	//  TPZFBMatrix  InnerProd(TPZFBMatrix &D );
	
	TPZFBMatrix &operator= (const TPZFBMatrix<TVar> & A );
	TPZFBMatrix operator+  (const TPZFBMatrix<TVar> & A ) const;
	TPZFBMatrix operator-  (const TPZFBMatrix<TVar> & A ) const;
	TPZFBMatrix &operator+=(const TPZFBMatrix<TVar> & A );
	TPZFBMatrix &operator-=(const TPZFBMatrix<TVar> & A );
	
	TPZFBMatrix operator*  (const TVar val ) const;
	TPZFBMatrix &operator*=(const TVar val );
	
	TPZFBMatrix operator-() const;
	
	long Dim() const     { return this->Rows(); }
	/** @brief Returns band size */
	long GetBand() const { return fBand; }
	/**
     * @brief Sets band size
     * @param newBand New band size
	 */
	int SetBand(const long newBand );
	
	/// Redimension the matrix preserving its elements
	int Resize(const long newRows,const long newCols );
	
	/// Redimension the matrix and make zero its elements
	int Redim(const long newRows,const long newCols );
	
	// Zeroes the elements of the matrix
	int Zero();
	
	void Transpose(TPZMatrix<TVar> *const T) const;
	int       Decompose_LU(std::list<long> &singular);
	int       Decompose_LU();
	
#ifdef OOPARLIB
	
	virtual long GetClassID() const   { return TFBMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	inline virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZFBMatrix" ); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const;
	
#endif
	
private:
	
	int Clear();
	
//	TPZVec<TVar> fElem;
	TVar *fElem;
	long  fBand;
};



/**************/
/*** PutVal ***/
template<class TVar>
inline int
TPZFBMatrix<TVar>::PutVal(const long row,const long col,const TVar& value )
{
	if ( (col+fBand >= row) && (col <= (row+fBand)) )
		fElem[((unsigned long) fBand * (2*row + 1)) + col ] = value;
	return( 1 );
}



/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &
TPZFBMatrix<TVar>::GetVal(const long row,const long col ) const {
	if ( (col+fBand >= row) && (col <= (row+fBand)) )
		return( fElem[ ((unsigned long)fBand * (2*row + 1)) + col ] );
	this->gZero = (TVar)0;
	return( this->gZero );
}

template<class TVar>
inline TVar &TPZFBMatrix<TVar>::operator()(const long row, const long col){  
	if( (col+fBand >= row) && (col <= (row+fBand)) )
	{
		return( fElem[ ((unsigned long)fBand * (2*row + 1)) + col ] );
	}
    DebugStop();
	this->gZero = (TVar)(0);
	return( this->gZero );
}
template<class TVar>
inline TVar &TPZFBMatrix<TVar>::s(const long row, const long col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}

#endif


