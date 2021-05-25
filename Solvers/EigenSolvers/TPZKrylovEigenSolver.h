#ifndef TPZKRYLOVEIGENSOLVER_H
#define TPZKRYLOVEIGENSOLVER_H
#include "TPZEigenSolver.h"


/** @brief Solvers for eigenvalue problems using Krylov methods.
    The eigenvalue problem is solved in the projected Krylov subspace 
    obtained by an Arnoldi iteration.*/
template<class TVar>
class TPZKrylovEigenSolver : public TPZEigenSolver<TVar>
{
public:
  using TPZEigenSolver<TVar>::TPZEigenSolver;

  TPZKrylovEigenSolver<TVar> * Clone() const override
  {return new TPZKrylovEigenSolver<TVar>(*this);}
  
  int SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  int SolveEigenProblem(TPZVec<CTVar> &w) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                   TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and does not calculates the eigenvectors
   * @param[out] w Stores the eigenvalues
   * @return Returns 1 if executed correctly
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w) override;
  
  //! Sets shift
  inline void SetShift(TVar s);
  //! Gets shift
  inline TVar Shift() const;
  /** @brief Sets number of Eigenpairs to calculate. 
      @note If the dimension of the Krylov subspace is insufficient, 
      it will be adjusted to `n`.
   */
  inline void SetNEigenpairs(int n);
  //! Gets number of Eigenpairs to calculate
  inline int NEigenpairs() const;
  //! Sets the dimension of the krylov subspace
  inline void SetKrylovDim(int n);
  //! Gets the dimension of the krylov subspace
  inline int KrylovDim() const;
  //! Sets tolerance for norm of the krylov vectors
  inline void SetTolerance(RTVar s);
  //! Gets tolerances
  inline RTVar Tolerance() const;

  inline void SetKrylovInitialVector(TPZFMatrix<TVar> q);

  inline TPZFMatrix<TVar> KrylovInitialVector() const;

  /** @brief Performs the Arnoldi iteration on a given matrix. 
      This iteration creates an orthonormal basis of the Krylov subspace of A. 
      The first vector of the Krylov basis can be set through 
      SetKrylovInitialVector method. The dimension of the Krylov subspace 
      can be set throught SetKrylovDim
      @param[in] A Matrix in which the algorithm is performed
      @param[out] Q Arnoldi vectors of A
      @param[out] H Representation of A in basis Q
      @return returns true if succeeded, false otherwise
  */
  [[nodiscard]] bool ArnoldiIteration(const TPZMatrix<TVar> &A,
                                      TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
                                      TPZFMatrix<TVar> &H);
protected:
  //! Shift to be applied
  TVar fShift{0};
  //! Number of Eigenpairs to calculate
  int fNEigenpairs{1};
  //! Dimension of the krylov subspace to calculate
  int fKrylovDim{1};
  //! Initial vector to be used to create Krylov subspace
  TPZFMatrix<TVar> fKrylovVector;
  //! Tolerance
  RTVar fTolerance{std::numeric_limits<RTVar>::epsilon()};
};


template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetShift(TVar s)
{
  fShift = s;
}
  
template<class TVar>
TVar TPZKrylovEigenSolver<TVar>::Shift() const
{
  return fShift;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetNEigenpairs(int n)
{
  if(n < 1) n = 1;
  fNEigenpairs = n;
  if(n > fKrylovDim){
    fKrylovDim = n;
    std::cout<< "Adjusted krylov dim to "<< fKrylovDim<<std::endl;
  }
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::NEigenpairs() const
{
  return fNEigenpairs;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetKrylovDim(int k)
{
  if(k<2){
    k = 2;
    std::cout<< "Adjusted krylov dim to "<< k<<std::endl;
  }
  fKrylovDim = k;
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::KrylovDim() const
{
  return fKrylovDim;
}


template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetTolerance(RTVar tol)
{
  fTolerance = tol;
}
  
template<class TVar>
RTVar TPZKrylovEigenSolver<TVar>::Tolerance() const
{
  return fTolerance;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetKrylovInitialVector(TPZFMatrix<TVar> q)
{
  fKrylovVector = q;
}

template<class TVar>
TPZFMatrix<TVar> TPZKrylovEigenSolver<TVar>::KrylovInitialVector() const
{
  return fKrylovVector;
}

extern template class TPZKrylovEigenSolver<float>;
extern template class TPZKrylovEigenSolver<double>;
extern template class TPZKrylovEigenSolver<std::complex<float>>;
extern template class TPZKrylovEigenSolver<std::complex<double>>;
#endif