//
// Created by Francisco Teixeira Orlandini on 18/05/2021
//

#ifndef TPZLAPACKEIGENSOLVER_H
#define TPZLAPACKEIGENSOLVER_H

#include "TPZEigenSolver.h"


template<class T>
class TPZSBMatrix;
/**
 * @ingroup solver
 * @brief  Defines an interface for using LAPACK to solve eigenvalue problems.
 * @note This class is only compatible with TPZFMatrix and TPZSBMatrix classes.
 */
template <typename TVar>
class TPZLapackEigenSolver : public TPZEigenSolver<TVar> {
  friend class TPZFMatrix<TVar>;
  friend class TPZSBMatrix<TVar>;
public:

  //!Default constructor.  
  TPZLapackEigenSolver();
  //!Copy constructor
  TPZLapackEigenSolver(const TPZLapackEigenSolver &copy) = default;
  //!Move constructor
  TPZLapackEigenSolver(TPZLapackEigenSolver &&rval) = default;
  //!Copy-assignment operator
  TPZLapackEigenSolver& operator=(const TPZLapackEigenSolver &copy) = default;
  //!Move-assignment operator
  TPZLapackEigenSolver& operator=(TPZLapackEigenSolver &&rval) = default;
  //!Destructor
  ~TPZLapackEigenSolver() = default;
  
  /** @name Eigen*/
  /** @{*/
  /**
   * @brief Solves the Ax=w*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns info param from LAPACK(0 if executed correctly)
   */
  int SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns info param from LAPACK(0 if executed correctly)
   */
  int SolveEigenProblem(TPZVec<CTVar> &w) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns info param from LAPACK(0 if executed correctly)
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and does not calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns info param from LAPACK(0 if executed correctly)
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w) override;
  /** 
      @brief Sets number of Eigenpairs to compute. 
      If not set, or set with n < 1, all pairs are returned.
      @note Since LAPACK only works will full matrices, this number
      will not affect computation, only the return value.
  */
  virtual void SetNEigenpairs(int n) override{
    this->fNEigenpairs = n;
  }
  /** @}*/
  //! Class identifier
  int ClassId() const override;
  //! Clone method
  TPZLapackEigenSolver<TVar>* Clone() const override;

  /**
     @brief Computes the eigenpairs of a Hessenberg matrix
     @param [in] A The Hessenberg matrix (it is only checked if Hessenberg in debug mode)
     @param[out] w The computed eigenvalues
     @param[out] eigenVectors The computed eigenvetors
     @return Returns info param from LAPACK(0 if executed correctly)
   */
  [[nodiscard]] int SolveHessenbergEigenProblem(TPZFMatrix<TVar> &A, TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors);
  /**
     @brief Computes the eigenvalues of a Hessenberg matrix
     @param [in] A The Hessenberg matrix (it is only checked if Hessenberg in debug mode)
     @param[out] w The computed eigenvalues
     @return Returns info param from LAPACK(0 if executed correctly)
   */
  [[nodiscard]] int SolveHessenbergEigenProblem(TPZFMatrix<TVar> &A, TPZVec<CTVar> &w);
protected:
  /*******************
   *    TPZFMATRIX    *
   *******************/
  /** @name EigenFMatrix */
  /** @{ */
  int SolveEigenProblem(TPZFMatrix<TVar> &A, TPZVec<CTVar> &w,
                        TPZFMatrix<CTVar> &eigenVectors, bool calcVectors=true);

  /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveEigenProblem(TPZFMatrix<TVar> &A, TPZVec<CTVar> &w);

  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param w Stores the eigenvalues
   * @param Stores the correspondent eigenvectors
   */
  int SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &A, TPZFMatrix<TVar> &B,
                                   TPZVec<CTVar> &w,
                                   TPZFMatrix<CTVar> &eigenVectors,
                                   bool calcVectors=true);
  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &A, TPZFMatrix<TVar> &B ,
                                   TPZVec<CTVar> &w);

  /** @} */
  /*******************
   *    TPSBMATRIX    *
   *******************/
  /** @name EigenSBMatrix */
  /** @{ */
  int SolveEigenProblem(TPZSBMatrix<TVar> &A, TPZVec<CTVar> &w,
                        TPZFMatrix<CTVar> &eigenVectors, bool calcVectors = true);

  /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveEigenProblem(TPZSBMatrix<TVar> &A, TPZVec<CTVar> &w);

  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param w Stores the eigenvalues
   * @param Stores the correspondent eigenvectors
   */
  int SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &A, TPZSBMatrix<TVar> &B ,
                                   TPZVec<CTVar> &w,
                                   TPZFMatrix<CTVar> &eigenVectors,
                                   bool calcVectors=true);
  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &A, TPZSBMatrix<TVar> &B ,
                                   TPZVec<CTVar> &w);
  /** @} */

  [[nodiscard]] int SolveHessenbergEigenProblem(TPZFMatrix<TVar> &A, TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors, bool calcEigenVectors);
};

extern template class TPZLapackEigenSolver<std::complex<float>>;
extern template class TPZLapackEigenSolver<std::complex<double>>;

extern template class TPZLapackEigenSolver<float>;
extern template class TPZLapackEigenSolver<double>;
#endif //TPZEIGENSOLVER_H
