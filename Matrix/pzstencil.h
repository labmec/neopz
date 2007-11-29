/******************************************************************************
 *
 * Class definition:    TPZStencilMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on sparse matrices stored by
 *                      stencils
 *
 * Solvers:             SOR
 *                      SSOR
 *
 *****************************************************************************/

#ifndef STENMATH
#define STENMATH

#include "pzmatrix.h"

class TPZFMatrix;

/**
Purpose:  Define operations on sparse matrices stored by stencils
*/
class TPZStencilMatrix : public TPZMatrix {

  public :

    TPZStencilMatrix( int rows, int cols );
  // sets up the StencilMatrix based on the stencil

  virtual ~TPZStencilMatrix();

  virtual int Rows() const;
  virtual int Cols() const;
  // Returns the rows or columns of this

  virtual const REAL &GetVal(const int row,const int col ) const;
  // Get the matrix entry at (row,col) without bound checking

  virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
		       const REAL alpha=1., const REAL beta = 0., const int opt = 0 , const int stride = 1) const;
  // computes z = beta * y + alpha * opt(this)*x
  //          z and x cannot overlap in memory

  virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
  // Print the matrix along with a identification title

  void SolveSOR( int &numiterations,const TPZFMatrix &rhs, TPZFMatrix &x,
		 TPZFMatrix *residual, TPZFMatrix &scratch,
		 const REAL overrelax, REAL &tol,
		 const int FromCurrent = 0,const int direction = 1 );

  void SetStencil( int stencilnumber, int inc, int *IA, REAL *A );
  // initiates Stencil number "stencilnumber" with the data

  void SetNodeStencils( int *stencilnumber );
  // associates the given stencil number with each row

 private:

  void IncreaseStencilPointers( int stencilnumber );

  void ComputeDiagonal();

  struct MPStencil {
    int fNumberOfItems;
    int fInc;
    int *fIA;
    REAL *fA;

    MPStencil( int inc, int *IA, REAL *A );
    ~MPStencil();

  } **fMystencils;

  int   fNumberOfStencilPointers;
  int   fRows, fCols;
  int  *fStencilNumbers;

  REAL *fDiag;
  int   fSolver;
  int   fMaxIterations;

  REAL fSORRelaxation;
};

inline int TPZStencilMatrix::Rows() const{
  return fRows;
}

inline int TPZStencilMatrix::Cols() const{
  return fCols;
}

#endif
