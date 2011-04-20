
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

class TPZFMatrix;
class TSFMatrix;



/**
@brief TPZSSpMatrix implement sparce symmetric matrix using a linked list of elements
*/
class TPZSSpMatrix : public TPZMatrix
{
 public:
  TPZSSpMatrix() : TPZMatrix(0,0)  {}
  TPZSSpMatrix(const int dim ) : TPZMatrix(dim,dim), fMat(dim, dim) {}
  TPZSSpMatrix(const TPZSSpMatrix & );
  
  CLONEDEF(TPZSSpMatrix)

  inline int    PutVal(const int row,const int col,const REAL&  element );
  inline const REAL & GetVal(const int row,const int col ) const;

  // Operadores com matrizes ESPARSAS simetricas.
  TPZSSpMatrix &operator= (const TPZSSpMatrix &A );
  TPZSSpMatrix operator+  (const TPZSSpMatrix &A ) const;
  TPZSSpMatrix operator-  (const TPZSSpMatrix &A ) const;
  TPZSSpMatrix &operator+=(const TPZSSpMatrix &A );
  TPZSSpMatrix &operator-=(const TPZSSpMatrix &A );

  // Operadores com matrizes GENERICAS.

  TPZSpMatrix operator+(const TPZSpMatrix &A ) const;
  TPZSpMatrix operator-(const TPZSpMatrix &A ) const;


  // Operadores com valores NUMERICOS.
  TPZSSpMatrix operator*  (const REAL v ) const;
  TPZSSpMatrix &operator*=(const REAL v );

  // Redimensiona a matriz, mas mantem seus elementos.
  int Resize(const int newDim ,const int )
    { fRow = fCol = newDim; return fMat.Resize( newDim, newDim ); }

  // Redimensiona a matriz e ZERA seus elementos.
  int Redim(const int newDim) { return Redim(newDim,newDim);}
  int Redim(const int newDim ,const int )
    { fRow = fCol = newDim; return fMat.Redim( newDim, newDim ); }

  // Zera a matriz
  int Zero()
    {return fMat.Zero();}


  /*** Resolucao de sistemas ***/

  int Decompose_Cholesky();  // Faz A = GGt.
  int Decompose_LDLt    ();  // Faz A = LDLt.
  int Decompose_Cholesky(std::list<int> &singular);  // Faz A = GGt.
  int Decompose_LDLt    (std::list<int> &singular);  // Faz A = LDLt.

  int Subst_Forward  ( TPZFMatrix *b ) const;
  //int Subst_Backward ( TPZMatrix *b );
  int Subst_LForward ( TPZFMatrix *b ) const;
  //int Subst_LBackward( TPZMatrix *b );
  int Subst_Diag     ( TPZFMatrix *b ) const;

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
  // Calcula o produto escalar entre as linhas 'row_i' e 'row_j'
  //  usando apenas os elementos pertencentes 'a colunas menores
  //  que 'k'. Ao retornar, as linhas 'row_i' e 'row_j' estarao
  //  posicionadas no elemento de coluna 'k' (ou onde ele deveria
  //  estar).
  REAL ProdEsc( TPZLink<TPZSpMatrix::TPZNode> *row_i,
		TPZLink<TPZSpMatrix::TPZNode> *row_j, int k );

//static int Error(const char *msg1,const char *msg2="" );


  TPZSpMatrix fMat;
};



/**************/
/*** PutVal ***/
//
//  Escreve um elemento na matriz, sem verificar fronteiras.
//  O valor da linha (row) deve ser maior ou igual ao da coluna
// (col).
//
inline int
TPZSSpMatrix::PutVal(const int r,const int c,const REAL&  value )
{
  // Inicializando row e col para trabalhar com a matriz
  //  triangular inferior.
  int row(r),col(c);
  if ( row < col )
    Swap( &row, &col );

  fDecomposed = 0;
  return( fMat.PutVal( row, col, value ) );
}



/**************/
/*** GetVal ***/
//
//  Le um elemento da matriz, sem verificar fronteiras. O valor
// da linha (row) deve ser maior ou igual ao da coluna (col).
//
inline const REAL &
TPZSSpMatrix::GetVal(const int r,const int c ) const
{
  // inicializando row e col para trabalhar com a matriz
  // triangular inferior.
  int row(r),col(c);
  if ( row < col )
    Swap( &row, &col );

  return( fMat.GetVal( row, col ) );
}

#endif
