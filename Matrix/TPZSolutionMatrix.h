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
  
  bool fIsComplex;//!< Whether it stores a complex or real solution
  TPZFMatrix<STATE> fRealMatrix;//!< Static storage for real matrix
  // TPZFMatrix<CSTATE> fComplexMatrix; //!< Static storage for complex matrix
  TPZBaseMatrix *fBaseMatrix;//!< Pointer for actual solution

public:
  /*! Constructor of empty TPZSolutionMatrix
      \param nrows - number of rows of the solution matrix
      \param ncols - number of cols of the solution matrix
      \param is_complex - whether the solution is complex or real
    */
  TPZSolutionMatrix(int nrows, int ncols, bool is_complex = false);

  //! Constructor taking a real matrix
  TPZSolutionMatrix(TPZFMatrix<STATE> &sol);
  // //! Constructor taking a complex matrix
  // TPZSolutionMatrix(TPZFMatrix<CSTATE> &sol);
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
  //! Get pointer to TPZBaseMatrix associated with the FEM solution
  inline TPZBaseMatrix *GetMatrixPtr() { return fBaseMatrix; }


  //the following doxygen directive will apply
  //this description for a group of functions
  
  //@{
  //! Get reference to real matrix (throws exception if solution is complex)
  inline const TPZFMatrix<STATE> &GetRealMatrix() const {
    if (fIsComplex)
      throw std::logic_error(
          "TPZSolutionMatrix has a complex matrix and GetRealMatrix was called");
    return fRealMatrix;
  }
  
  inline TPZFMatrix<STATE> &GetRealMatrix() {
    if (fIsComplex)
      throw std::logic_error(
          "TPZSolutionMatrix has a complex matrix and GetRealMatrix was called");
    return fRealMatrix;
  }
  //@}

  // const TPZFMatrix<CSTATE> &GetComplexMatrix() const { return fRealMatrix;
  // }
  // TPZFMatrix<CSTATE> &GetComplexMatrix() { return fRealMatrix; }

  //! ClassId method
  int ClassId() const override{
    return Hash("TPZSolutionMatrix");
  }
  //! Read method
  void Read(TPZStream &buf, void *context) override;
  //! Write method (only the current matrix is written)
  void Write(TPZStream &buf, int withclassid) const override;
};

#endif