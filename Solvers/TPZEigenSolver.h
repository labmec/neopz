//
// Created by Francisco Teixeira Orlandini on 18/05/2021
//

#ifndef TPZEIGENSOLVER_H
#define TPZEIGENSOLVER_H

#include "TPZSolver.h"

#include "pzfmatrix.h"
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
  /** @}*/
  //! Class identifier
  int ClassId() const override;
  //! Resets Matrices
  void ResetMatrix() override;
protected:
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

extern template class TPZEigenSolver<float>;
extern template class TPZEigenSolver<double>;
extern template class TPZEigenSolver<long double>;

extern template class TPZEigenSolver<std::complex<float>>;
extern template class TPZEigenSolver<std::complex<double>>;
extern template class TPZEigenSolver<std::complex<long double>>;
#endif //TPZEIGENSOLVER_H
