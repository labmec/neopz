//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatred.h
//
// Class:  TPZMatRed
//
// Obs.:   Subestruturacao simples de um sistema de equacoes.
//
//				[K00][U0] + [K01][U1] = [F0]
//				[K10][U0] + [K11][U1] = [F1]
//
// Versao: 04 / 1996.
//


#ifndef _TMATREDHH_
#define _TMATREDHH_

//#include "tintvec.h"



#include "pzmatrix.h"
#include "pzreal.h"

#ifdef OOPARLIB
#include "pzsaveable.h"
#include "pzmatdefs.h"
#endif


class TPZFMatrix;


/**
 * Implements a simple substructuring of a linear system of equations
 *
 *			[K00][U0] + [K01][U1] = [F0]
 *			[K10][U0] + [K11][U1] = [F1]
 *@ingroup matrix
*/
class TPZMatRed :public TPZMatrix {
  //	friend class TPZMatrix;

 public:
  /**
   *Simple constructor
   */
  TPZMatRed ();


  TPZMatRed(const  int dim,const int dim00  ) ;
  /**
   *Simple destructor
   */
  ~TPZMatRed ();

  virtual int IsSimetric() const;
  // returns 1 or 0 depending on whether the fK00 matrix is zero or not

  // Put and Get values without bounds checking
  // these methods are faster than "Put" e "Get" if DEBUG is defined
  int    PutVal(const int row,const int col,const REAL& value );
  const REAL &GetVal(const int row,const int col ) const;

  /**
   * This method will zero all submatrices associated with this reducable matrix class
   **/
  virtual int Zero();

  /**
   * Sets the matrix pointer of the upper left matrix to K00
   * @param k00 pointer to an upper left matrix
   */
  void SetK00(TPZMatrix *const k00);
  /**
   * Copies the F vector in the internal data structure
   * @param F vector containing data to stored in current object
   */
  void SetF(const TPZFMatrix & F);
  void SetDecomposeType(const DecomposeType t = ELU){fDecomposeType=t;}


  /**
   * Computes the reduced version of the right hand side
   * [F1]=[F1]-[K10][A00^-1][F0]
   */
  const TPZFMatrix & F1Red();

  /**
   * Computes the K11 reduced
   * [K11]=[K11]-[K10][A00^-1][A01]
   */
  const TPZFMatrix & K11Red();

  /**
   * Returns the second vector, inverting K11
   * @param F contains second vector
   */
  void U1(TPZFMatrix & F);

  /**
   * Computes the complete vector based on the solution u1.
   * @param U1 right hand side
   * @param result contains the result of the operation
   */
  void UGlobal(const TPZFMatrix & U1, TPZFMatrix & result);
  //  TPZFMatrix U0(TPZMatrix *u1 = NULL);



  void Print(const char *name = NULL, std::ostream &out = std::cout,const MatrixOutputFormat = EFormatted ) const;

  int Substitution(TPZFMatrix *right_side) const;

  /**
	*Redim: Set the dimension of the complete matrix
	*and reduced matrix
	**/
	int Redim(int dim, int dim00); //Cesar 19/12/00

 private:

//static int Error(const char *msg ,const char *msg2 = "");

  /**
   * If fK00 is simetric, only part of the matrix is accessible to external
   * objects. Simetrizes copies the data of the matrix to make its data simetric
   */
  void Simetrize();

  /**
   * Swaps the row and column index
   * @param row Row number
   * @param col Column number
   */
  static void Swap(int *row, int *col);

  /**
   * Stiffnes matrix
   */
  TPZMatrix *fK00;
  /**
   * Full Stiffnes matrix
   */
  TPZFMatrix *fK01,
    *fK10,
    *fK11;
  /**
   * Right hand side or force matrix
   */
  TPZFMatrix *fF0,
    *fF1;

  DecomposeType  fDecomposeType;

  /**
   * Stores matricess fKij dimensions
   */
  int fDim0,fDim1;

  /**
   * Is true if [K00^-1][F0] has been calculated and overwritten on [F0]
   */
  char fF0IsComputed;

  /**
   * Is true if [K00^-1][K01] has been calculated and overwritten on [K01]
   */
  char fK01IsComputed;

  /**
   * fK11IsReduced is true if [K11]=[K11]-[K10][A00^-1][A01] exists
   */
  char fK11IsReduced;

  /**
   * fF1IsReduced is true if  [F1]=[F1]-[K10][A00^-1][F0] exists
   */
  char fF1IsReduced;
};

/************/
/*** Swap ***/
/* Modificacao por Philippe Devloo (insercao de inline )*/
inline void
TPZMatRed::Swap( int *a, int *b )
{
  int aux = *a;
  *a = *b;
  *b = aux;
}


#endif
