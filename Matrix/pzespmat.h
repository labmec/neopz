
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tespmat.hh
//
// Class:  TPZSpMatrix
//
// Obs.:   Implementa matrizes esparsas. A implementacao
//         e' feita atraves de listas ligadas.
//
// Versao: 04 / 1996.
//


#ifndef TESPMATH
#define TESPMATH



#include "pzlink.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef WORKPOOL

#include "pzworkpool.h"

#endif

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

class TPZFMatrix;



/**
 * Defines sparce matrix class. \n
 * Stores data in a sparce storage scheme
 * @ingroup matrix
 */
class TPZSpMatrix : public TPZMatrix
{
  friend class TPZSSpMatrix;

 public:
  
  /**
     @struct TPZNode
     Defines a node
     @param elem Element on node
     @param col Column number
  */
   struct TPZNode
  {
    REAL elem;
    int    col;
  };

 public:
  /**
   *Simple constructor
   */
  TPZSpMatrix ()
    : TPZMatrix( 0, 0 )  { fElem = NULL; }
  /**
     Constructor with initialization parameters
     @param rows Number of rows
     @param cols Number of columns
  */
  TPZSpMatrix ( int rows, int cols );
  /**
     Copy constructor
     @param A Model object
  */
  TPZSpMatrix (const TPZSpMatrix &A )
    : TPZMatrix( A.Rows(), A.Cols() )  { fCopy( &A ); }
  /**
     Simple destructor
  */
  ~TPZSpMatrix();

  int    Put(const int row,const int col,const REAL & value );
  const REAL &Get(const int row,const int col ) const;

  // Nao verifica limites da matriz (e' mais rapido).
  int    PutVal(const int row,const int col,const REAL & element );
  const REAL &GetVal(const int row,const int col ) const;

  // Operadores com matrizes ESPARSAS NAO simetricas.
  //@{
  /**
     Generic overloaded operator
  */ 
  TPZSpMatrix &operator= (const TPZSpMatrix &A );
  TPZSpMatrix operator+  (const TPZSpMatrix &A ) const;
  TPZSpMatrix operator-  (const TPZSpMatrix &A ) const;
  //  TPZSpMatrix operator*  ( TPZSpMatrix &A );
  TPZSpMatrix &operator+=(const TPZSpMatrix &A );
  TPZSpMatrix &operator-=(const TPZSpMatrix &A );

  // Operadores com matrizes GENERICAS.
  TPZSpMatrix &operator=(const TPZMatrix &A );
  //@}
  virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
		       const REAL alpha ,const REAL beta = 0.,const int opt = 0,const int stride = 1) const;


  // Operadores com valores NUMERICOS.
  //@{
  /**
     Numerical values operator
  */
  TPZSpMatrix operator*  (const REAL v ) const;
  TPZSpMatrix &operator*=(const REAL v );

  TPZSpMatrix operator-() const { return operator*(-1.0); }
  //@}
  // Desaloca todos elementos da matriz, ou seja, ZERA a matriz.
  TPZSpMatrix &Reset();

  // Redimensiona a matriz, mas mantem seus elementos.
  int Resize(const int newRows,const int newCols );
  /**
     Redimensions current matrix keeping its elements
     @param newDim New matrix dimensio
  */
  int Resize(const int newDim )   { return Resize( newDim, newDim ); }

  // Redimensiona a matriz e ZERA seus elementos.
  int Redim(const int newRows,const int newCols );
  /**
     Redimensions matrix deleting all its elements
     @param newDim New matrix dimension
  */
  int Redim(const int newDim )    { return Redim( newDim, newDim ); }

  // Zera os elementos da matriz
  int Zero();

  /*** Resolucao de sistemas ***/

#ifdef OOPARLIB

  virtual long GetClassID() const   { return TSPMATRIX_ID; }
  virtual int Unpack( TReceiveStorage *buf );
  static TSaveable *Restore(TReceiveStorage *buf);
  inline virtual int Pack( TSendStorage *buf ) const;
  virtual char *ClassName() const   { return( "TPZSpMatrix" ); }
  virtual int DerivedFrom(const long Classid) const;
  virtual int DerivedFrom(const char *classname) const;

#endif

 protected:
  //  int fAdd(const REAL &val );// operator +=
  /**
   *  @see TPZSpMatrix::operator+=
   *  @link TPZSpMatrix::operator+=
   */
  int fAdd(const TPZSpMatrix *const A );               // operator +=
  int fAdd(const TPZSpMatrix *A, TPZSpMatrix *B ); // operator +

  int fSub(const TPZSpMatrix *const A );               // operator -=
  int fSub( TPZSpMatrix *A, TPZSpMatrix *B ); // operator -

  //  int fCopy( REAL val );                // operator =
  int fCopy(const TPZSpMatrix *const A );              // operator =

  int fMult( REAL val );                // operator *


  // Troca o valor das duas variaveis.
  inline void Swap( int *a, int *b );

  int  Clear()               { delete [] fElem; return ( 1 );}

  // Calcula o produto escalar entre as linhas 'row_i' e 'row_j'
  //  usando apenas os elementos pertencentes 'a colunas menores
  //  que 'k'. Ao retornar, as linhas 'row_i' e 'row_j' estarao
  //  posicionadas no elemento de coluna 'k' (ou onde ele deveria
  //  estar).
  /**
     Computes dot product with respect to lines row_i e row_j using only elements
     that belongs to columns less than k. Afterwards lines row_i and row_j would 
     be positioned on the element of column k. (where it should be).
     @param row_i Ith row to be used
     @param row_j Jth row to be used
  */
  REAL ProdEsc( TPZLink<TPZNode> *row_i, TPZLink<TPZNode> *row_j,
		int k );


  TPZLink<TPZNode> *fElem;
#ifdef WORKPOOL
  TPZWorkPool fWp;
#endif

 private:

  int    Error (const char *msg1,const char *msg2="" ) const;
};



/*** Swap ***/

inline void
TPZSpMatrix::Swap( int *a, int *b )
{
  int c = *a;
  *a = *b;
  *b = c;
}

#endif
