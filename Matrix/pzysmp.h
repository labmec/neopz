/******************************************************************************
 *
 * Class definition:    TPZFYsmpMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on general sparse matrices stored
 *                      in the (old) Yale Sparse Matrix Package format.
 *
 * Solvers:             SOR
 *
 *****************************************************************************/

#ifndef YSMPMATH
#define YSMPMATH

#include "pzmatrix.h"

class TPZFMatrix;

/**
Define operations on general sparse matrices stored
in the (old) Yale Sparse Matrix Package format.
*/
class TPZFYsmpMatrix : public TPZMatrix {

  public :

    TPZFYsmpMatrix(const int rows,const int cols );
  // sets up the StencilMatrix based on the stencil

  virtual ~TPZFYsmpMatrix();


  virtual const REAL &GetVal(const int row,const int col ) const;
  // Get the matrix entry at (row,col) without bound checking


  void PutVal(const int row, const int col, double Value);

  virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
		       const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 )const;
  // computes z = beta * y + alpha * opt(this)*x
  //          z and x cannot overlap in memory

  virtual void SetData( int *IA, int *JA, REAL *A );
  // Pass the data to the class.

  virtual void Print(const char *title, ostream &out = cout , const MatrixOutputFormat form = EFormatted) const;
  // Print the matrix along with a identification title

  void SolveSOR(int &numiterations, const TPZFMatrix &rhs, TPZFMatrix &x,
		TPZFMatrix *residual, TPZFMatrix &scratch,
		const REAL overrelax, REAL &tol,
		const int FromCurrent = 0,const int direction = 1 ) const;    

  /**
       * Add a contribution of a stiffness matrix
  	* putting it on destination indexes position
       */
  virtual void AddKelOld(
  		TPZFMatrix & elmat //! Member stiffness matrix beeing added
  		, TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
  		);    

  void AddKel(TPZFMatrix & elmat, TPZVec<int> & destinationindex);    

  void Multiply(TPZFYsmpMatrix & B, TPZFYsmpMatrix & Res);    

  virtual int Zero();


 private:

  void ComputeDiagonal();
protected:
  int  *fIA;
  int  *fJA;
  REAL *fA;

  REAL *fDiag;

  int   fSymmetric;

  //    int   fSolver;
  //    int   fMaxIterations;

  //    REAL  fSORRelaxation;
protected:

  /**
   * Implements a initialization method for the sparse structure. It sets the initial value for the fIA and fJA.
   * 
   * -fIA will contain the initial positions for all the equations
   * -fJA will contain (-1) on all its positions
   * 
   * -fA will contain 0 on all its value 
   */
  void InitializeData();
};



inline void TPZFYsmpMatrix::SetData( int *IA, int *JA, REAL *A ) {
  // Pass the data to the class.
  fIA = IA;
  fJA = JA;
  fA  =  A;
  ComputeDiagonal();
}

#endif
