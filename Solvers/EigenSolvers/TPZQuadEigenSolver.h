#ifndef TPZQUADEVPSOLVER_H
#define TPZQUADEVPSOLVER_H
#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolverBase.h"

/**
   @brief This class implements a solver for the generalised quadratic EVP
   -b^2 Mu + i b Lu + Ku = 0

   where M,L and K are matrices, i is the imaginary unit,
   b is the eigenvalue and u the associated eigenvector.

   For now it only uses the shift and invert strategy for computing
   eigenvalues around a given value.

   The quadratic problem is thus reduced to a generalised EVP

   Q (u)  = ib P (u)
     (w)         (w)

   where

   w = i b u

   and

   Q = (-K 0)
       ( 0 I)

   and
   
   P = (L M)
       (I 0)
       
   This is a Krylov based solver that will use Arnoldi iteration.

   For reference, see

   Waveguide Propagation Modes and Quadratic Eigenvalue Problems

   by
   
   André Nicolet and Christophe Geuzaine
 */
template<class TVar>
class TPZQuadEigenSolver : public TPZEigenSolver<TVar>,
                           public TPZKrylovEigenSolverBase<TVar>{
public:
    //!Default constructor.  
    TPZQuadEigenSolver() = default;
    //!Copy constructor
    TPZQuadEigenSolver(const TPZQuadEigenSolver &copy) = default;
    //!Move constructor
    TPZQuadEigenSolver(TPZQuadEigenSolver &&rval) = default;
    //!Copy-assignment operator
    TPZQuadEigenSolver& operator=(const TPZQuadEigenSolver &copy) = default;
    //!Move-assignment operator
    TPZQuadEigenSolver& operator=(TPZQuadEigenSolver &&rval) = default;
    //!Destructor
    virtual ~TPZQuadEigenSolver() = default;

    inline int Solve(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override
    {
        return this->SolveImpl(w,eigenVectors,true);
    }
    inline int Solve(TPZVec<CTVar> &w) override
    {
        TPZFMatrix<CTVar> ev;
        return this->SolveImpl(w,ev,false);
    }

    void SetMatrixK(TPZAutoPointer<TPZMatrix<TVar>> K);
    void SetMatrixL(TPZAutoPointer<TPZMatrix<TVar>> L);
    void SetMatrixM(TPZAutoPointer<TPZMatrix<TVar>> M);

    TPZAutoPointer<TPZMatrix<TVar>> MatrixK() {return fK;}
    TPZAutoPointer<TPZMatrix<TVar>> MatrixL() {return fL;}
    TPZAutoPointer<TPZMatrix<TVar>> MatrixM() {return fM;}
    
    void SetMatrices(TPZAutoPointer<TPZMatrix<TVar>> K,
                     TPZAutoPointer<TPZMatrix<TVar>> L,
                     TPZAutoPointer<TPZMatrix<TVar>> M);


    void SetTarget(TVar s) override {
        TPZEigenSolver<TVar>::SetTarget(s);
        SetShift(s);
    }
    //! Sets shift of the shift and invert spectral transform
    void SetShift(TVar s){fShift = s;}
    //! Gets shift of the shift and invert spectral transform
    TVar Shift() const {return fShift;}    
    //! Applies (maybe matrix-free) operator on a given vector
    void ApplyOperator(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &res) const override;
    //! Number of rows of the eigenvector
    [[nodiscard]] int64_t NRows() const override{
        return this->fM->Rows();
    }
    //! Algebraic size (twice the number of rows)
    [[nodiscard]] int64_t SystemSize() const override{
        return this->NRows()*2;
    }

    TPZQuadEigenSolver<TVar> * Clone() const override;
    void ResetMatrix() override;
    void SetNEigenpairs(int n) override;
protected:
    //!Computes the only matrix to be decomposed delta = -(s^2 M + sL + K), where s is the shift
    void ComputeAndDecomposeDelta();

    void PreSolve() override;
    void TransformEigenvalues(TPZVec<CTVar> &w) override;
    //! Matrices of the EVP
    TPZAutoPointer<TPZMatrix<TVar>> fM,fL,fK;
    //! Only matrix to be decomposed, delta = -(s^2 M + sL + K), where s is the shift
    TPZAutoPointer<TPZMatrix<TVar>> fDelta;
    //! Shift for the shift invert transform
    TVar fShift;
    //! Temporary vector for avoiding repeated memory allocations
    mutable TPZFMatrix<TVar> fScratch;
};


extern template class TPZQuadEigenSolver<float>;
extern template class TPZQuadEigenSolver<double>;

extern template class TPZQuadEigenSolver<std::complex<float>>;
extern template class TPZQuadEigenSolver<std::complex<double>>;

#endif