#ifndef TPZKRYLOVEIGENSOLVERBASE_H
#define TPZKRYLOVEIGENSOLVERBASE_H
#include <pzreal.h>

#include <pzfmatrix.h>

template<class T>
class TPZAutoPointer;

/** @brief Base class for solvers for eigenvalue problems using Krylov methods.
    It basically implements a Krylov based solver using the Arnoldi iteration.
    See TPZKrylovEigenSolver and TPZQuadEigenSolver for usage*/
template<class TVar>
class TPZKrylovEigenSolverBase
{
public:
  TPZKrylovEigenSolverBase() = default;
  //we need a virtual dtor
  ~TPZKrylovEigenSolverBase() = default;
  //thus we need copy/move constructors
  TPZKrylovEigenSolverBase(const TPZKrylovEigenSolverBase<TVar> &cp) = default;
  TPZKrylovEigenSolverBase& operator=(const TPZKrylovEigenSolverBase<TVar> &cp) = default;
  TPZKrylovEigenSolverBase(TPZKrylovEigenSolverBase<TVar> &&rval) = default;
  TPZKrylovEigenSolverBase& operator=(TPZKrylovEigenSolverBase<TVar> &&rval) = default;
  /** @brief Sets the dimension of the Krylov subspace. 
      @note If not set, defaults to `10*nev`, where `nev` is the number 
      of sought eigenvalues.*/
  inline void SetKrylovDim(int d);
  //! Gets the dimension of the Krylov subspace
  inline int KrylovDim() const;
  //! Sets tolerance for norm of the Krylov vectors
  inline void SetTolerance(RTVar s);
  //! Gets tolerances
  inline RTVar Tolerance() const;
  /** @brief Sets the first vector to be used in the Krylov subspace
      @note The input vector will be normalized. If not set, a random 
      vector will be generated.
  */
  inline void SetKrylovInitialVector(TPZFMatrix<TVar> q);
  //! Gets the first vector used in the Krylov subspace
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
  [[nodiscard]] bool ArnoldiIteration(TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
                                      TPZFMatrix<TVar> &H);

  //! Applies (maybe matrix-free) operator on a given vector
  virtual void ApplyOperator(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &res) const = 0;
  //! Algebraic system size (normally corresponds to number of rows)
  [[nodiscard]] virtual int64_t SystemSize() const {return this->NRows();}
  //! Number of rows of the eigenvectors
  [[nodiscard]] virtual int64_t NRows() const = 0;
protected:
  //! Dimension of the Krylov subspace to calculate
  int64_t fKrylovDim{-1};
  //! Initial vector to be used to create Krylov subspace
  TPZFMatrix<TVar> fKrylovVector;
  //! Tolerance
  RTVar fTolerance{std::numeric_limits<RTVar>::epsilon()};
  //! Implementation of Solve methods. Need not be implemented in derived classes.
  virtual int SolveImpl(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors, bool computeVectors);

  //! Derived classes can implement operations to be performed before solving
  virtual void PreSolve() {};
  //! Derived classes should implement how to transform eigenvalues (i.e., based on spectral transformations)
  virtual void TransformEigenvalues(TPZVec<CTVar>&w) = 0;
};


template<class TVar>
void TPZKrylovEigenSolverBase<TVar>::SetKrylovDim(int k)
{
  fKrylovDim = k;
}

template<class TVar>
int TPZKrylovEigenSolverBase<TVar>::KrylovDim() const
{
  return fKrylovDim;
}


template<class TVar>
void TPZKrylovEigenSolverBase<TVar>::SetTolerance(RTVar tol)
{
  fTolerance = tol;
}
  
template<class TVar>
RTVar TPZKrylovEigenSolverBase<TVar>::Tolerance() const
{
  return fTolerance;
}

template<class TVar>
void TPZKrylovEigenSolverBase<TVar>::SetKrylovInitialVector(TPZFMatrix<TVar> q)
{
  fKrylovVector = q;
}

template<class TVar>
TPZFMatrix<TVar> TPZKrylovEigenSolverBase<TVar>::KrylovInitialVector() const
{
  return fKrylovVector;
}
#endif
