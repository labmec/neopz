/**
 * @file
 * @brief Contains TPZSSpMatrix class which implements sparce symmetric matrix using a linked list of elements.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tsespmat.hh
//
// Class:  TPZSSpMatrix
//
// Deriva: TSimMatrix
//
// Obs.:   Implementa matrizes esparsas simetricas. A implementacao
//         e' feita atraves de listas ligadas.
//
// Versao: 04 / 1996.
//


#ifndef TSESPMATH
#define TSESPMATH

#include "pzlink.h"
#include "pzespmat.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzsfulmat.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TSFMatrix;

/**
 * @brief Implements sparce symmetric matrix using a linked list of elements. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSSpMatrix : public TPZMatrix<TVar>
{
public:
	TPZSSpMatrix() : TPZMatrix<TVar>(0,0)  {}
	TPZSSpMatrix(const int dim ) : TPZMatrix<TVar>(dim,dim), fMat(dim, dim) {}
	TPZSSpMatrix(const TPZSSpMatrix<TVar> & );
	
	CLONEDEF(TPZSSpMatrix)
	
	inline int    PutVal(const int row,const int col,const TVar&  element );
	inline const TVar & GetVal(const int row,const int col ) const;
	
	/// Operators with SYMMETRIC sparse matrices.
	// @{
	TPZSSpMatrix &operator= (const TPZSSpMatrix<TVar> &A );
	TPZSSpMatrix operator+  (const TPZSSpMatrix<TVar> &A ) const;
	TPZSSpMatrix operator-  (const TPZSSpMatrix<TVar> &A ) const;
	TPZSSpMatrix &operator+=(const TPZSSpMatrix<TVar> &A );
	TPZSSpMatrix &operator-=(const TPZSSpMatrix<TVar> &A );
	// @}
	
	/// Operators with Generic matrices.
	// @{
	TPZSpMatrix<TVar> operator+(const TPZSpMatrix<TVar> &A ) const;
	TPZSpMatrix<TVar> operator-(const TPZSpMatrix<TVar> &A ) const;
	// @}
	
	/// Operators with numeric values.
	// @{
	TPZSSpMatrix operator*  (const TVar v ) const;
	TPZSSpMatrix &operator*=(const TVar v );
	// @}
	
	/// Resize the array but keeps its entirety.
	int Resize(const int newDim ,const int )
    { this->fRow = this->fCol = newDim; return fMat.Resize( newDim, newDim ); }
	
	/// Resize the array and resets its entirety.
	int Redim(const int newDim) { return Redim(newDim,newDim);}
	int Redim(const int newDim ,const int )
    { this->fRow = this->fCol = newDim; return fMat.Redim( newDim, newDim ); }
	
	// Zeroes all the elements
	int Zero()
    {return fMat.Zero();}
	
	
	/*** @brief To solve linear systems ***/
	// @{
	int Decompose_Cholesky();  // Faz A = GGt.
	int Decompose_LDLt    ();  // Faz A = LDLt.
	int Decompose_Cholesky(std::list<int> &singular);  // Faz A = GGt.
	int Decompose_LDLt    (std::list<int> &singular);  // Faz A = LDLt.
	
	int Subst_Forward  ( TPZFMatrix<TVar> *b ) const;
	//int Subst_Backward ( TPZMatrix<>*b );
	int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
	//int Subst_LBackward( TPZMatrix<>*b );
	int Subst_Diag     ( TPZFMatrix<TVar> *b ) const;
	// @}
	
#ifdef OOPARLIB
	
	virtual long GetClassID() const        { return TSSPMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual char *ClassName() const   { return( "TPZSSpMatrix"); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const; // a class with name classname
	
#endif
	
private:
	/** @brief Calcula o produto escalar entre as linhas 'row_i' e 'row_j'
	 *  usando apenas os elementos pertencentes 'a colunas menores
	 *  que 'k'.
	 *
	 * Ao retornar, as linhas 'row_i' e 'row_j' estarao
	 *  posicionadas no elemento de coluna 'k' (ou onde ele deveria
	 *  estar).
	 */
	TVar ProdEsc( TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_i,
				 TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_j, int k );
	
	TPZSpMatrix<TVar> fMat;
};


/**************/
/*** PutVal ***/
//
//  Escreve um elemento na matriz, sem verificar fronteiras.
//  O valor da linha (row) deve ser maior ou igual ao da coluna
// (col).
//
template<class TVar>
inline int
TPZSSpMatrix<TVar>::PutVal(const int r,const int c,const TVar&  value )
{
	// Inicializando row e col para trabalhar com a matriz
	//  triangular inferior.
	int row(r),col(c);
	if ( row < col )
		this->Swap( &row, &col );
	
	this->fDecomposed = 0;
	return( fMat.PutVal( row, col, value ) );
}

/**************/
/*** GetVal ***/
//
//  Le um elemento da matriz, sem verificar fronteiras. O valor
// da linha (row) deve ser maior ou igual ao da coluna (col).
//
template<class TVar>
inline const TVar &
TPZSSpMatrix<TVar>::GetVal(const int r,const int c ) const
{
	// inicializando row e col para trabalhar com a matriz
	// triangular inferior.
	int row(r),col(c);
	if ( row < col )
		this->Swap( &row, &col );
	
	return( fMat.GetVal( row, col ) );
}

#endif
