//
// Created by Francisco Teixeira Orlandini on 18/05/2021
//

#ifndef TPZEIGENSOLVER_H
#define TPZEIGENSOLVER_H

#include "TPZSolver.h"
#include "pzfmatrix.h"
#include "TPZEigenSort.h"


template <class T>
class TPZKrylovEigenSolverBase;


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
  virtual int Solve(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) = 0;

  /**
   * @brief Solves the EVP (according to IsGeneralised()) 
   and does not calculate the eigenvectors.
   The default implementation relies on calls to SolveGeneralisedEigenProblem
   and SolveEigenProblem, to be implemented in base classes.
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  virtual int Solve(TPZVec<CTVar> &w) = 0;
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
  /** @brief Decides criterium for sorting the obtained eigenvalues. 
      @note By default it is set to TPZEigenSort::EAbsAscending .*/
  inline void SetEigenSorting(TPZEigenSort ord);
  //! Returns criterium for sorting the obtained eigenvalues
  inline TPZEigenSort EigenSorting() const;
  /** @}*/
  //! Class identifier
  int ClassId() const override;

  friend class TPZKrylovEigenSolverBase<TVar>;
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
  /**
   * @brief Stores the computed eigenvalues
   */
  TPZVec<CTVar> fEigenvalues;
  /**
   * @brief Stores the computed eigenvectors
   */
  TPZFMatrix<CTVar> fEigenvectors;
  //! Sorting order of the eigenvalues
  TPZEigenSort fEigenSort{TPZEigenSort::AbsAscending};
  //! Target eigenvalue
  TVar fTarget{0};
};

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
