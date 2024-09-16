#ifndef TPZKRYLOVEIGENSOLVER_H
#define TPZKRYLOVEIGENSOLVER_H
#include "TPZKrylovEigenSolverBase.h"
#include "TPZLinearEigenSolver.h"
#include "TPZSpectralTransform.h"


/** @brief Solvers for linear eigenvalue problems using Krylov methods.
    The eigenvalue problem is solved in the projected Krylov subspace 
    obtained by an Arnoldi iteration. 
    See TPZSpectralTransform for possible spectral transformations*/
template<class TVar>
class TPZKrylovEigenSolver : public TPZLinearEigenSolver<TVar>,
  public TPZKrylovEigenSolverBase<TVar>
{
public:
  using TPZLinearEigenSolver<TVar>::TPZLinearEigenSolver;

  TPZKrylovEigenSolver<TVar> * Clone() const override
  {return new TPZKrylovEigenSolver<TVar>(*this);}

  /** @name BasicUsage */
  /** @{*/
  //! Set target eigenvalue and adjust spectral transform if needed
  void SetTarget(TVar target) override;
  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  int SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @return Returns 1 if executed correctly
   */
  int SolveEigenProblem(TPZVec<CTVar> &w) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                   TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and does not calculates the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @return Returns 1 if executed correctly
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w) override;

  //! Sets spectral transformation to be used (it creates an internal copy)
  void SetSpectralTransform(TPZSpectralTransform<TVar> &s);
  //! Gets spectral transformation to be used
  inline TPZAutoPointer<TPZSpectralTransform<TVar>> SpectralTransform();
  /** 
      @brief Sets number of Eigenpairs to calculate. 
      @note If the dimension of the Krylov subspace is insufficient, 
      it will be adjusted to `n`.
  */
  inline void SetNEigenpairs(int n) override;
  /** @}*/


  //! Applies (maybe matrix-free) operator on a given vector
  void ApplyOperator(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &res) const override;
  //! Algebraic size (number of rows)
  [[nodiscard]] int64_t SystemSize() const override {return NRows();}
  //! Number of rows of the eigenvector
  [[nodiscard]] int64_t NRows() const override;
protected:
  //! Operations to be performed at the beginning of SolveImpl call
  void PreSolve() override;
  //! Transforms eigenvalues based on spectral transforms
  void  TransformEigenvalues(TPZVec<CTVar>&w) override;
  //! Spectral Transformation
  TPZAutoPointer<TPZSpectralTransform<TVar>> fST;
  //! Set target of spectral transformation equal to user-defined target
  void AdjustTargetST();
  //! whether the target has been set by user
  bool fUserTarget{false};
};


template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetSpectralTransform(TPZSpectralTransform<TVar> &s)
{
  fST = s.Clone();
}

template<class TVar>
TPZAutoPointer<TPZSpectralTransform<TVar> >
TPZKrylovEigenSolver<TVar>::SpectralTransform()
 {
  return fST;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetNEigenpairs(int n)
{
  if(n < 1) n = 1;
  this->fNEigenpairs = n;
  if(n > this->fKrylovDim){
    this->fKrylovDim = 10 * n;
    std::cout<< "Adjusted Krylov dim to "<<this->fKrylovDim<<std::endl;
  }
}

extern template class TPZKrylovEigenSolver<float>;
extern template class TPZKrylovEigenSolver<double>;
extern template class TPZKrylovEigenSolver<std::complex<float>>;
extern template class TPZKrylovEigenSolver<std::complex<double>>;
#endif
