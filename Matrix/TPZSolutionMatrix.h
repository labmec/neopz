/**
 * @file
 * @brief Contains declaration of the TPZSolutionMatrix class, which provides an abstraction for dealing with real/complex FEM solution matrix.
 */

#ifndef TPZSOLUTIONMATRIX_H
#define TPZSOLUTIONMATRIX_H

#include "pzreal.h"//STATE
#include "pzfmatrix.h"//TPZFMatrix<T>,TPZBaseMatrix
#include "Hash/TPZHash.h"//Hash

class TPZSolutionMatrix : public TPZSavable{
private:

  //! Enum for solution type
  enum ESolType{ EReal=1,EComplex=2,EUndefined=3};
  ESolType fSolType;//!< Type of solution
  TPZFMatrix<STATE> fRealMatrix;//!< Static storage for real matrix
  // TPZFMatrix<CSTATE> fComplexMatrix; //!< Static storage for complex matrix
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
  // //! Constructor taking a complex matrix
  // TPZSolutionMatrix(const TPZFMatrix<CSTATE> &sol);
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
  //! Conversion function to TPZBaseMatrix&
  operator TPZBaseMatrix& ();
  
  // //@{
  // //! Conversion function to TPZFMatrix<CSTATE>&. Throws error if incompatible
  // operator TPZFMatrix<CSTATE>& ();
  // operator const TPZFMatrix<CSTATE>& () const;
  // //@}
  
  //! Number of Rows of the solution
  inline int64_t Rows() const { return fBaseMatrix->Rows(); }
  //! Number of cols of the solution
  inline int64_t Cols() const { return fBaseMatrix->Cols(); }
  //! Redim the solution \ref matrix "Matrix"
  inline int Redim(const int64_t r, const int64_t c) {
    return fBaseMatrix->Redim(r, c);
  }
  //! Resize the solution \ref matrix "Matrix"
  inline int Resize(const int64_t r, const int64_t c) {
    return fBaseMatrix->Resize(r, c);
  }
  //! Prints the matrix in an std::ostream
  void Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted){
    return fBaseMatrix->Print(name,out,EFormatted);
  }
  //! Zeroes the matrix
  int Zero() {
    return fBaseMatrix->Zero();
  }
  //! ClassId method
  int ClassId() const override{
    return Hash("TPZSolutionMatrix");
  }
  //! Read method
  void Read(TPZStream &buf, void *context) override;
  //! Write method (only the current matrix is written)
  void Write(TPZStream &buf, int withclassid) const override;
};


extern template
TPZSolutionMatrix &TPZSolutionMatrix::operator=<STATE>(const TPZFMatrix<STATE> &mat);
// extern template
// TPZSolutionMatrix &TPZSolutionMatrix::operator=<CSTATE>(const TPZFMatrix<CSTATE> &mat);
#endif