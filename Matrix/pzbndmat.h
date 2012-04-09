/**
 * @file
 * @brief Contains TPZFBMatrix class which defines a non symmetric banded matrix.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tbndmat.hh
//
// Class:  TPZFBMatrix
//
// Obs.:   Implementa matrizes cheias (normais).
//
// Versao: 12 / 1994.
//

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
 */
template<class TVar>
class TPZFBMatrix : public TPZMatrix<TVar>
{
	
public:
	virtual int Substitution(TPZFMatrix<TVar> *B) const;
	/**
	 * @brief Simple constructor
	 */
	TPZFBMatrix ();
	/**
     @brief Constructor
     @param dim Initial dimension of banded matrix
     @param band_width Initial band width of banded matrix
	 */
	TPZFBMatrix (const int dim,const int band_width = 0 );
	/**
     @brief Copy constructor
	 */
	TPZFBMatrix (const TPZFBMatrix<TVar> & );
	
	CLONEDEF(TPZFBMatrix)
	/**
     @brief Simple destructor
	 */
	~TPZFBMatrix();
	
	int    Put(const int row,const int col,const TVar& value );
	const TVar &Get(const int row,const int col ) const;
	
	TVar &operator()(const int row, const int col);
	virtual TVar &s(const int row, const int col);
	//estos metodos nao verificam a existencia do elemento
	//sao mas rapidos que Put e Get
	inline int    PutVal(const int row,const int col,const TVar& value );
	inline const TVar &GetVal(const int row,const int col ) const;
	
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
	
	int Dim() const     { return this->Rows(); }
	/**
     @brief Returns band size
	 */
	int GetBand() const { return fBand; }
	/**
     @brief Sets band size
     @param newBand New band size
	 */
	int SetBand(const int newBand );
	
	// Redimensiona a matriz, mas mantem seus elementos.
	// Nao muda o tamanho da banda!
	/// Redimension the matrix preserving its elements
	int Resize(const int newRows,const int newCols );
	
	// Redimensiona a matriz e ZERA seus elementos.
	// Nao muda o tamanho da banda!
	/// Redimension the matrix and make zero its elements
	int Redim(const int newRows,const int newCols );
	
	// Zera os elementos da matriz
	int Zero();
	
	void Transpose(TPZMatrix<TVar> *const T) const;
	int       Decompose_LU(std::list<int> &singular);
	int       Decompose_LU();
	
#ifdef OOPARLIB
	
	virtual long GetClassID() const   { return TFBMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	inline virtual int Pack( TSendStorage *buf ) const;
	virtual char *ClassName() const   { return( "TPZFBMatrix" ); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const;
	
#endif
	
private:
	
	
	//static  int Error(const char *msg1,const char *msg2="" );
	int Clear();
	
	TVar *fElem;
	int  fBand;
};



/**************/
/*** PutVal ***/
template<class TVar>
inline int
TPZFBMatrix<TVar>::PutVal(const int row,const int col,const TVar& value )
{
	if ( (col+fBand >= row) && (col <= (row+fBand)) )
		fElem[ fBand * (2*row + 1) + col ] = value;
	return( 1 );
}



/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &
TPZFBMatrix<TVar>::GetVal(const int row,const int col ) const {
	if ( (col+fBand >= row) && (col <= (row+fBand)) )
		return( fElem[ fBand * (2*row + 1) + col ] );
	this->gZero = 0.;
	return( this->gZero );
}

template<class TVar>
inline TVar &TPZFBMatrix<TVar>::operator()(const int row, const int col){  
	if( (col+fBand >= row) && (col <= (row+fBand)) )
	{
		return( fElem[ fBand * (2*row + 1) + col ] );
	}
    DebugStop();
	this->gZero = 0.;
	return( this->gZero );
}
template<class TVar>
inline TVar &TPZFBMatrix<TVar>::s(const int row, const int col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}

#endif


