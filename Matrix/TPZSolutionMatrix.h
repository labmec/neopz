/**
 * @file
 * @brief Contains declaration of the TPZSolutionMatrix class, which provides an abstraction for dealing with real/complex FEM solution matrix.
 */

#ifndef TPZSOLUTIONMATRIX_H
#define TPZSOLUTIONMATRIX_H

#include "pzreal.h"//STATE
#include "pzfmatrix.h"//TPZFMatrix<T>,TPZBaseMatrix
#include "Hash/TPZHash.h"//Hash


//definition for debugging access to fBaseMatrix
//#define PZ_SOLMAT_DEBUG

class TPZSolutionMatrix : public TPZSavable{
private:

  //! Enum for solution type
  ESolType fSolType;//!< Type of solution
  TPZFMatrix<STATE> fRealMatrix;//!< Static storage for real matrix
  TPZFMatrix<CSTATE> fComplexMatrix; //!< Static storage for complex matrix
  TPZBaseMatrix *fBaseMatrix;//!< Pointer for actual solution

public:
  //! Default constructor
  TPZSolutionMatrix();
  //! Constructor defining matrix type
  TPZSolutionMatrix(bool is_complex);
  /*! Constructor of empty TPZSolutionMatrix
      \param nrows - number of rows of the solution matrix
      \param ncols - number of cols of the solution matrix
      \param is_complex - whether the solution is complex or real
    */
  TPZSolutionMatrix(int nrows, int ncols, bool is_complex = false);

  //! Constructor taking a real matrix
  TPZSolutionMatrix(const TPZFMatrix<STATE> &sol);
  //! Constructor taking a complex matrix
  TPZSolutionMatrix(const TPZFMatrix<CSTATE> &sol);
  //! Copy constructor
  TPZSolutionMatrix(const TPZSolutionMatrix &);
  //! Move constructor (deleted)
  TPZSolutionMatrix(TPZSolutionMatrix &&) = delete;
  //! Destructor
  ~TPZSolutionMatrix() = default;
  //! Copy assignment operator
  TPZSolutionMatrix &operator=(const TPZSolutionMatrix &);
  //! Move assignment operator
  TPZSolutionMatrix &operator=(TPZSolutionMatrix &&) = delete;
  //! Assignment operator from TPZFMatrix<T>. Throws error if incompatible
  template<class TVar>
  TPZSolutionMatrix &operator=(const TPZFMatrix<TVar> &mat);

  //@{
  //! Conversion function to TPZFMatrix<STATE>&. Throws error if incompatible
  operator TPZFMatrix<STATE>& ();
  operator const TPZFMatrix<STATE>& () const;
  //@}

  //@{
  //! Conversion function to TPZFMatrix<CSTATE>&. Throws error if incompatible
  operator TPZFMatrix<CSTATE>& ();
  operator const TPZFMatrix<CSTATE>& () const;
  //@}
  
  //! Conversion function to TPZBaseMatrix&
  operator TPZBaseMatrix& ();


  //@{
  //!Arithmetic operators. Throws error if incompatible.
  template<class TVar>
  TPZSolutionMatrix &operator+=(const TPZFMatrix<TVar> &mat);

  TPZSolutionMatrix &operator+=(const TPZSolutionMatrix &sol);
  //@}
  
  /*! Sets a solution with nrows where the first sol.Rows() are the solution from sol.
    Since the TPZAnalysis only stores the solution associated with the independent
    equations, in TPZCompMesh we still need extra rows for the solution associated
    with the dependent equations. This method provided a way for doing this operation
    with only one memory allocation.*/
  void ExpandAndSetSol(const TPZSolutionMatrix & sol, const int64_t nrows);
  
  //! Number of Rows of the solution
  inline int64_t Rows() const {
#ifdef PZ_SOLMAT_DEBUG
    if(!fBaseMatrix) DebugStop();
#endif
    return fBaseMatrix->Rows();
  }
  //! Number of cols of the solution
  inline int64_t Cols() const {
#ifdef PZ_SOLMAT_DEBUG
    if(!fBaseMatrix) DebugStop();
#endif
    return fBaseMatrix->Cols();
  }
  //! Redim the solution \ref matrix "Matrix"
  inline int Redim(const int64_t r, const int64_t c) {
#ifdef PZ_SOLMAT_DEBUG
    if(!fBaseMatrix) DebugStop();
#endif
    if(fBaseMatrix) return fBaseMatrix->Redim(r, c);
    return 0;
  }
  //! Resize the solution \ref matrix "Matrix"
  inline int Resize(const int64_t r, const int64_t c) {
#ifdef PZ_SOLMAT_DEBUG
    if(!fBaseMatrix) DebugStop();
#endif
    if(fBaseMatrix) return fBaseMatrix->Resize(r, c);
    return 0;
  }
  //! Prints the matrix in an std::ostream
  inline void Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted){
#ifdef PZ_SOLMAT_DEBUG
    if(!fBaseMatrix) DebugStop();
#endif
    if(fBaseMatrix) return fBaseMatrix->Print(name,out,EFormatted);
  }
  //! Zeroes the matrix
  inline int Zero() {
#ifdef PZ_SOLMAT_DEBUG
    if(!fBaseMatrix) DebugStop();
#endif
    if(fBaseMatrix) return fBaseMatrix->Zero();
    return 0;
  }
  //! ClassId method
  int ClassId() const override{
    return Hash("TPZSolutionMatrix");
  }
  //! Read method
  void Read(TPZStream &buf, void *context) override;
  //! Write method (only the current matrix is written)
  void Write(TPZStream &buf, int withclassid) const override;
  friend STATE Norm(const TPZSolutionMatrix &A);
};

inline STATE Norm(const TPZSolutionMatrix &A) {
  if (A.fSolType == EReal)
    return Norm(A.fRealMatrix);
  else if (A.fSolType == EComplex){
    return Norm(A.fComplexMatrix);
  }
  
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"Type is not set. Aborting...\n";
  DebugStop();
  return -1;
}

extern template
TPZSolutionMatrix &TPZSolutionMatrix::operator=<STATE>(const TPZFMatrix<STATE> &mat);
extern template
TPZSolutionMatrix &TPZSolutionMatrix::operator+=<STATE>(const TPZFMatrix<STATE> &mat);
extern template
TPZSolutionMatrix &TPZSolutionMatrix::operator=<CSTATE>(const TPZFMatrix<CSTATE> &mat);
extern template
TPZSolutionMatrix &TPZSolutionMatrix::operator+=<CSTATE>(const TPZFMatrix<CSTATE> &mat);
#endif