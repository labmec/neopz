//
// Created by Francisco Teixeira Orlandini on 18/05/2021
//

#ifndef TPZEIGENSOLVER_H
#define TPZEIGENSOLVER_H

#include "TPZSolver.h"
#include "pzfmatrix.h"

//! Sorting method for calculated eigenvalues
enum class TPZEigenSort{
  AbsAscending,/*!< Ascending magnitude*/
  AbsDescending,/*!< Descending magnitude*/
  RealAscending,/*!< Ascending real part*/
  RealDescending,/*!< Descending real part*/
  ImagAscending,/*!< Ascending imaginary part*/
  ImagDescending,/*!< Descending imaginary part*/
  TargetRealPart,/*!< Real part closest to target*/
  TargetImagPart,/*!< Imaginary part closest to target*/
  TargetMagnitude/*!< Magnitude closest to target*/
};

/**
* @ingroup solver
* @brief  Defines an interface for  eigenvalue problems solvers.
*/
template <typename TVar>
class TPZEigenSolver : public TPZSolver {
public:
  //!Default constructor.  
  TPZEigenSolver() = default;
  //!Copy constructor
  TPZEigenSolver(const TPZEigenSolver &copy) = default;
  //!Move constructor
  TPZEigenSolver(TPZEigenSolver &&rval) = default;
  //!Copy-assignment operator
  TPZEigenSolver& operator=(const TPZEigenSolver &copy) = default;
  //!Move-assignment operator
  TPZEigenSolver& operator=(TPZEigenSolver &&rval) = default;
  //!Destructor
  virtual ~TPZEigenSolver() = default;

  /** @name Eigen*/
  /** @{*/
  /**
   * @brief Solves the EVP (according to IsGeneralised()) and 
   calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  inline int Solve(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors);

  /**
   * @brief Solves the EVP (according to IsGeneralised()) 
   and does not calculate the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  inline int Solve(TPZVec<CTVar> &w);
  /**
   * @brief Solves the Ax=w*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  virtual int SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) = 0;

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  virtual int SolveEigenProblem(TPZVec<CTVar> &w) = 0;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  virtual int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                           TPZFMatrix<CTVar> &eigenVectors) = 0;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and does not calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  virtual int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w) = 0;
  /** 
      @brief Sets number of Eigenpairs to compute. 
      @note In some solvers, n<1 is equal to all.
  */
  virtual void SetNEigenpairs(int n) = 0;
  //! Gets number of Eigenpairs to calculate
  inline int NEigenpairs() const{
    return fNEigenpairs;
  }
  /** @brief Sets target eigenvalue
      @note If there is a shift of origin or shift and invert 
      spectral transformation associated with
      the solver, then setting a target will also set the
      corresponding shift
  */
  inline virtual void SetTarget(TVar t);
  
  inline TVar Target() const;
  //!Gets the Matrix A
  inline TPZAutoPointer<TPZMatrix<TVar>> MatrixA(){
    return fMatrixA;
  }
  //!Gets the Matrix B(for generalised eigenvalue problems)
  inline TPZAutoPointer<TPZMatrix<TVar>> MatrixB(){
    return fMatrixB;
  }
  //!Sets the Matrix A
  inline void SetMatrixA(TPZAutoPointer<TPZMatrix<TVar>> mat){
    fMatrixA = mat;
  }
  //!Sets the Matrix B (for generalised eigenvalue problems)
  inline void SetMatrixB(TPZAutoPointer<TPZMatrix<TVar>> mat){
    fMatrixB = mat;
  }
  //!Whether the solver is set for a generalised eigenvalue problem
  inline bool IsGeneralised() const{
    return fIsGeneralised;
  }
  //!Configure the solver to solve a generalised eigenvalue problem
  virtual void SetAsGeneralised(bool isGeneralised);
  /** @brief Decides criterium for sorting the obtained eigenvalues. 
      @note By default it is set to TPZEigenSort::EAbsAscending .*/
  inline void SetEigenSorting(TPZEigenSort ord);
  //! Returns criterium for sorting the obtained eigenvalues
  inline TPZEigenSort EigenSorting() const;
  /** @}*/
  //! Class identifier
  int ClassId() const override;
  //! Resets Matrices
  void ResetMatrix() override;
protected:
  /**
     @brief Sort the calculated eigenvalues and return a vector with
     size given by the NEigenvalues() method
     @param [in/out] eigenvalues
     @param [out] vector of indices relating the unsorted and sorted eigenvalues.
  */
  void SortEigenvalues(TPZVec<CTVar> &w, TPZVec<int> &indices);
  //! Number of Eigenpairs to calculate
  int fNEigenpairs{-1};
  /** @brief Whether to solve the eigenvalue problem
   *   is generalised (Ax=uBx) or not (Ax=ux)*/
  bool fIsGeneralised{false};
  /**
   * @brief Stores the computed eigenvalues
   */
  TPZVec<CTVar> fEigenvalues;
  /**
   * @brief Stores the computed eigenvectors
   */
  TPZFMatrix<CTVar> fEigenvectors;

  /** @brief Container classes */
  TPZAutoPointer<TPZMatrix<TVar>> fMatrixA{nullptr};

  /** @brief Container classes */
  TPZAutoPointer<TPZMatrix<TVar>> fMatrixB{nullptr};
  //! Sorting order of the eigenvalues
  TPZEigenSort fEigenSort{TPZEigenSort::AbsAscending};
  //! Target eigenvalue
  TVar fTarget{0};
};

template<class TVar>
int TPZEigenSolver<TVar>::Solve(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors)
{
  if(IsGeneralised()) return SolveGeneralisedEigenProblem(w,eigenVectors);
  else return SolveEigenProblem(w,eigenVectors);
}

template<class TVar>
int TPZEigenSolver<TVar>::Solve(TPZVec<CTVar> &w)
{
  if(IsGeneralised()) return SolveGeneralisedEigenProblem(w);
  else return SolveEigenProblem(w);
}

template<class TVar>
TPZEigenSort TPZEigenSolver<TVar>::EigenSorting() const
{
  return fEigenSort;
}

template<class TVar>
void TPZEigenSolver<TVar>::SetEigenSorting(TPZEigenSort ord)
{
  fEigenSort = ord;
}

template<class TVar>
void TPZEigenSolver<TVar>::SetTarget(TVar target)
{
  fTarget = target;
}

template<class TVar>
TVar TPZEigenSolver<TVar>::Target() const
{
  return fTarget;
}

extern template class TPZEigenSolver<float>;
extern template class TPZEigenSolver<double>;
extern template class TPZEigenSolver<long double>;

extern template class TPZEigenSolver<std::complex<float>>;
extern template class TPZEigenSolver<std::complex<double>>;
extern template class TPZEigenSolver<std::complex<long double>>;
#endif //TPZEIGENSOLVER_H
