
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tfullmat.hh
//
// Class:  TPZFMatrix
//
// Obs.:   Implements matrix classes (normais).
//
// Versao: 04 / 1996.
//
#ifndef _TMATRIXHH_
#include "pzmatrix.h"
#endif


#ifndef _TFULLMATRIXH_
#define _TFULLMATRIXH_


#include <iostream>

using namespace std;

#include <math.h>

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

/******************************/

class TPZTempFMatrix;
template <class T>
class TPZVec;
/**
 * Non abstract class which defines full matrices
 * @ingroup matrix
 */
class TPZFMatrix : public TPZMatrix {

 public:
  /**
   *Simple constructor
   */
  TPZFMatrix () : TPZMatrix( 0, 0 ), fGiven(0) { fElem = NULL; fSize = 0;}
  /**
     Constructor with initialization parameters
     @param rows Initial number of rows
     @param cols Number of columns
     @param buf Who knows !!!
     @param size Who knows too !!!
  */
  TPZFMatrix (const int rows ,const int columns = 1, REAL * buf = NULL,const int size = 0 );
  /**
     Constructor with initialization parameters
     @param rows Initial number of rows
     @param cols Number of columns
     @param val Inital value fill all elements
  */
  TPZFMatrix (const int rows ,const int columns,const REAL & val );
  //@{
  /**
     Copy constructor
     @param refmat Used as a model for current object
  */
  TPZFMatrix (const TPZFMatrix & );
  TPZFMatrix(const TPZMatrix &refmat); // copy the elements one by one
  //@}
  /**
     Constructor that uses a temporary matrix
  */
  TPZFMatrix(TPZTempFMatrix );
  /**
     Simple destructor
  */
  virtual  ~TPZFMatrix();

  int PutVal(const int row,const int col,const REAL & value );
  const REAL &GetVal(const int row,const int col ) const;
  
  virtual REAL &s(const int row, const int col);

  REAL &g(const int row, const int col) const;
  /**
   * Performs a right hand side assemblage
   * @param rhs Load vector
   * @param destinantion Destine index on current matrix
   */
  void AddFel(TPZFMatrix &rhs,TPZVec<int> &destination);
  /**
   * Performs a right hand side assemblage
   * @param rhs Load vector
   * @param source Source index on rhs
   * @param destinantion Destine index on current matrix
   */
  void AddFel(TPZFMatrix &rhs,TPZVec<int> &source, TPZVec<int> &destination);
  virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
		       const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 ) const;
 
  //@{
  /**
     Generic operator with REAL type
  */
  REAL &operator()(const int row,const int col);
  REAL &operator()(const int row);
  //@}
  /**
     @name FULL
     Operations with FULL matrices
  */
  //@{
  /**
     Generic operator with FULL matrices
  */
  TPZFMatrix &operator= (const TPZFMatrix &A );
  TPZFMatrix &operator= ( TPZTempFMatrix A);
  TPZTempFMatrix operator+  (const TPZFMatrix &A ) const;
  //TPZTempFMatrix operator+ (TPZTempFMatrix A);
  TPZTempFMatrix operator-  (const TPZFMatrix &A ) const;
  TPZTempFMatrix operator*  (const TPZFMatrix &A ) const;
  TPZFMatrix &operator+=(const TPZFMatrix &A );
  //    TPZFMatrix &operator+=(TPZTempFMatrix A );
  TPZFMatrix &operator-=(const TPZFMatrix &A );
  //@}
  /**
   * Performs an ZAXPY operation being *this += alpha * p
   * @param alpha Being alpha on above opereation
   * @param p Being p on above operation
   */
  void ZAXPY(const REAL alpha, const TPZFMatrix &p);
  /**
   * Performs an operation *this = this * beta + z
   * @param beta Being beta on above opereation
   * @param z Being z on above operation
   */
  void TimesBetaPlusZ(const REAL beta, const TPZFMatrix &z);

  // Operations with matrices GENERICAS.
  /**
     @name Generics
     Generic operators with matrices
  */
  //@{
  /**
     Generic operator with matrices
  */
  TPZFMatrix &operator= (const TPZMatrix &A );
  TPZTempFMatrix operator+  (const TPZMatrix &A ) const;
  TPZTempFMatrix operator-  (const TPZMatrix &A ) const;
  TPZTempFMatrix operator*  (const TPZMatrix &A ) const;
  TPZFMatrix &operator+=(const TPZMatrix &A );
  TPZFMatrix &operator-=(const TPZMatrix &A );
  //@}
  // Operations with values NUMERICOS.
  
  /**
     @name Numerics
     Numeric operations with matrices
  */
  //@{
  /**
     Numeric operator with matrices
  */
  TPZFMatrix &operator= (const REAL val );
  TPZTempFMatrix operator+  (const REAL val ) const;
  TPZTempFMatrix operator-  (const REAL val ) const;// { return operator+( -val ); }
  TPZTempFMatrix operator*  (const REAL val ) const;
  TPZFMatrix &operator+=(const REAL val );
  TPZFMatrix &operator-=(const REAL val )  { return operator+=( -val ); }
  TPZFMatrix &operator*=(const REAL val );

  TPZTempFMatrix operator-() const;// { return operator*( -1.0 ); }
  //@}
  //  void Input( istream & in = cin );

  // Redimension a matrix, but maintain your elements.
  int Resize(const int newRows,const int wCols );

  // Redimension a matrix and ZERO your elements.
  int Redim(const int newRows,const int newCols );

  // Zero os elements
  int Zero();

  void Transpose(TPZMatrix *const T) const;
  /**
     @see TPZMatrix::Transpose
  */
  void Transpose();

  /*** Solve some systems ***/

  int Decompose_LU();
  int Substitution( TPZFMatrix *B ) const;


  //routines to send and receive messages

#ifdef OOPARLIB
  virtual long GetClassID() const   { return TFMATRIX_ID; }
  virtual int Unpack( TReceiveStorage *buf );
  static TSaveable *Restore(TReceiveStorage *buf);
  inline virtual int Pack( TSendStorage *buf ) const;
  virtual char *ClassName() const   { return( "TPZFMatrix" ); }
  virtual int DerivedFrom(const long Classid) const; // returns  true if the object
  //  belongs to a class which is derived from a class
  //  with id classid
  virtual int DerivedFrom(const char *classname) const; // a class with name classname
#endif

 private:

  int Error(const char *msg1,const char *msg2=0 ) const;
  int Clear();

  REALPtr fElem;
  REALPtr fGiven;
  int fSize;
};


REAL Dot(const TPZFMatrix &A,const TPZFMatrix &B);

REAL Norm(const TPZFMatrix &A);

/**************/
/*** PutVal ***/
inline int TPZFMatrix::PutVal(const int row, const int col,const REAL & value ) {
  fElem[ ((unsigned)col) * Rows() + row ] = value;
  return( 1 );
}



/**************/
/*** GetVal ***/
inline const REAL &TPZFMatrix::GetVal( const int row, const int col ) const {
  return( fElem[ ((unsigned)col) * Rows() + row ] );
}


inline REAL &TPZFMatrix::operator()( const int row, const int col) {
#ifndef NODEBUG
  if(row >= Rows() || row<0 || col >= Cols() || col<0) {
    Error("TPZFMatrix::operator()","Index out of bounds");
    return gZero;
  }
#endif
  return *(fElem+col*fRow+row);
}

inline REAL &TPZFMatrix::s(const int row, const int col) {
  // verificando se o elemento a inserir esta dentro da matriz
  return operator()(row,col);
}


inline REAL &TPZFMatrix::g( const int row, const int col) const {
#ifdef DEBUG
  if(row >= Rows() || row<0 || col >= Cols() || col<0) {
    Error("TPZFMatrix::operator()","Index out of bounds");
    return gZero;
  }
#endif
  return *(fElem+col*fRow+row);
}

inline REAL &TPZFMatrix::operator()(const int row) {
#ifdef DEBUG
  if(row >= Rows() || row<0) {
    Error("TPZFMatrix::operator()","Index out of bounds");
    return gZero;
  }
#endif
  return *(fElem+row);
}

/**************************/
/*** Operations Global ***/



//inline TPZFMatrix &TPZFMatrix::operator+=(TPZTempFMatrix A) {
//	return (*this) += A.Object();
//}

inline REAL Norm(const TPZFMatrix &A) {
  return sqrt(Dot(A,A));
}


#endif


