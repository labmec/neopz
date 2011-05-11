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

#define TPZMATRED_VERYSPARSE_ID 28291103
#define TPZMATRED_FMATRIX_ID 28291102

#include "pzmatrix.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzsolve.h"
#include "tpzverysparsematrix.h"

#ifdef OOPARLIB
#include "pzsaveable.h"
#include "pzmatdefs.h"
#endif

class TPZFMatrix;

/**
 * @brief Implements a simple substructuring of a linear system of equations

 *			[K00][U0] + [K01][U1] = [F0]
 *			[K10][U0] + [K11][U1] = [F1]
 *@brief Implements a matrix composed of 4 submatrices
 *@ingroup matrix
 */
/// Implements a matrix composed of 4 submatrices
template<class TSideMatrix = TPZFMatrix>
class TPZMatRed: public TPZMatrix
{
    // @ingroup matrix
public:
	
	friend class TPZMatRed<TPZFMatrix>;
	friend class TPZMatRed<TPZVerySparseMatrix>;
  /**
   *Simple constructor
   */
  TPZMatRed();

  /**
   * Constructor with 2 parameters
   * dim assumes the value of n1+n2
   * dim00 equals n1
   */
  TPZMatRed(const int dim, const int dim00);

	template<class TSideCopy>
	TPZMatRed<TSideMatrix>(const TPZMatRed<TSideCopy> &cp): TPZMatrix(cp), fK11(cp.fK11), fK01(cp.fK01), fK10(cp.fK10), fF0(cp.fF0), fF1(cp.fF1)
	{
		fDim0=cp.fDim0;
		fDim1=cp.fDim1;
		fF0IsComputed=cp.fF0IsComputed;
		fK11IsReduced=cp.fK11IsReduced;
		fK01IsComputed = cp.fK01IsComputed;
		fF1IsReduced=cp.fF1IsReduced;
		fIsReduced = cp.fIsReduced;
		fSolver = cp.fSolver;
		
		if(cp.fK00) fK00 = cp.fK00;
	}
	
	

  CLONEDEF(TPZMatRed)
  /**
   *Simple destructor
   */
  ~TPZMatRed();

  /**
   * returns 1 or 0 depending on whether the fK00 matrix is zero or not
   */
  virtual int IsSimetric() const;

  /**
   * changes the declared dimension of the matrix to fDim1
   */
  void SetReduced()
  {
    TPZMatrix::Resize(fDim1, fDim1);
	  fIsReduced = 1;
  }

  /**
   * Put and Get values without bounds checking
   * these methods are faster than "Put" e "Get" if DEBUG is defined
   */
  virtual int PutVal(const int row, const int col, const REAL& value);
  virtual const REAL &GetVal(const int row, const int col) const;
  virtual REAL &s(int row, int col);

  /**
   * This method will zero all submatrices associated with this reducable matrix class
   **/
  virtual int Zero();

  /**
   * Sets the matrix pointer of the upper left matrix to K00
   * @param k00 pointer to an upper left matrix
   */
  void SetK00(TPZAutoPointer<TPZMatrix> K00);

  TPZAutoPointer<TPZMatrix> K00()
  {
    return fK00;
  }
  TSideMatrix &K01()
  {
    return fK01;
  }
  TSideMatrix &K10()
  {
    return fK10;
  }
  void SetSolver(TPZAutoPointer<TPZMatrixSolver> solver);
  /**
   * Copies the F vector in the internal data structure
   * @param F vector containing data to stored in current object
   */
  void SetF(const TPZFMatrix & F);
	
	/**
	 * Indicate whether F0 needs to be reduced or not
	 */
	void SetF0IsComputed(bool value)
	{
		fF0IsComputed = value;
	}
	/*
	 * Indicate that the value of F1 has been reduced
	 */
	void SetF1IsReduced(bool value)
	{
		fF1IsReduced = value;
	}

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
	void UGlobal2(TPZFMatrix & U1, TPZFMatrix & result);
  //  TPZFMatrix U0(TPZMatrix *u1 = NULL);


  /**
   * Prints the object data structure
   */
  void Print(const char *name = NULL, std::ostream &out = std::cout,
      const MatrixOutputFormat = EFormatted) const;


  /**
   * Redim: Set the dimension of the complete matrix
   * and reduced matrix
   */
  int Redim(int dim, int dim00); //Cesar 19/12/00


  /**
   * It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
   * @param x Is x on the above operation
   * @param y Is y on the above operation
   * @param z Is z on the above operation
   * @param alpha Is alpha on the above operation
   * @param beta Is beta on the above operation
   * @param opt Indicates if is Transpose or not
   * @param stride Indicates n/N where n is dimension of the right hand side
   * vector and N is matrix dimension
   */
  void MultAdd(const TPZFMatrix &x, const TPZFMatrix &y, TPZFMatrix &z,
      const REAL alpha, const REAL beta, const int opt, const int stride) const;

  /**
   * If fK00 is simetric, only part of the matrix is accessible to external
   * objects. Simetrizes copies the data of the matrix to make its data simetric
   */
  void Simetrize();

  /**
   * template class TPZMatRed<TPZVerySparseMatrix>;
template class TPZMatRed<TPZFMatrix>;
   *
   */
  /**
   * Saveable methods
   */
  int ClassId() const;

  //TPZMATRED_FMATRIX_ID
  virtual void Write(TPZStream &buf, int withclassid);
  virtual void Read(TPZStream &buf, void *context);

private:

  //static int Error(const char *msg ,const char *msg2 = "");

  /**
   * Swaps the row and column index
   * @param row Row number
   * @param col Column number
   */
  static void Swap(int *row, int *col);

  /**
   * Stiffnes matrix
   */
  TPZAutoPointer<TPZMatrix> fK00;

  /// Solution method for inverting fK00
  TPZAutoPointer<TPZMatrixSolver> fSolver;
  /**
   * Full Stiffnes matrix
   */
  TSideMatrix fK11;
  TSideMatrix fK01, fK10;

  /**
   * Right hand side or force matrix
   */
  TPZFMatrix fF0, fF1;


  /**
   * Stores matricess fKij dimensions
   */
  int fDim0, fDim1;

  /**
   * Is true if the declared dimension of the matrix is fDim0
   */
  char fIsReduced;

  /**
   * Is true if [K00^-1][F0] has been calculated and overwritten on [F0]
   */
  char fF0IsComputed;

  /**
   * Is true if [K00^-1][KO1] has been computed and overwritten [K01]
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
template<class TSideMatrix>
inline void TPZMatRed<TSideMatrix>::Swap(int *a, int *b)
{
  int aux = *a;
  *a = *b;
  *b = aux;
}


#endif
