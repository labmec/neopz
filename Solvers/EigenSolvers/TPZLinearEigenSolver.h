#ifndef TPZLINEAREIGENSOLVER_H
#define TPZLINEAREIGENSOLVER_H
#include "TPZEigenSolver.h"


template<class T>
class TPZPardisoSolver;

/**
   @brief This class implements a solver for the
   linear (possibily generalised) Eigenvalue Problem
   Au = l u (Au = l Bu)

   where A (and B) is a matrix, l the eigenvalue and u the
   associated eigenvector.
*/
template<class TVar>
class TPZLinearEigenSolver : public TPZEigenSolver<TVar>{
public:
    //!Default constructor.  
    TPZLinearEigenSolver() = default;
    //!Copy constructor
    TPZLinearEigenSolver(const TPZLinearEigenSolver &copy) = default;
    //!Move constructor
    TPZLinearEigenSolver(TPZLinearEigenSolver &&rval) = default;
    //!Copy-assignment operator
    TPZLinearEigenSolver& operator=(const TPZLinearEigenSolver &copy) = default;
    //!Move-assignment operator
    TPZLinearEigenSolver& operator=(TPZLinearEigenSolver &&rval) = default;
    //!Destructor
    virtual ~TPZLinearEigenSolver() = default;

    /** @name Eigen*/
    /** @{*/
    /**
     * @brief Solves the EVP (according to IsGeneralised()) and 
     calculates the eigenvectors
     * @param[out] w Stores the eigenvalues
     * @return Returns 1 if executed correctly
     */
    inline int Solve(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override{
        if(this->IsGeneralised()) return SolveGeneralisedEigenProblem(w, eigenVectors);
        else return SolveEigenProblem(w, eigenVectors);
    }

    /**
     * @brief Solves the EVP (according to IsGeneralised()) 
     and does not calculate the eigenvectors.
     The default implementation relies on calls to SolveGeneralisedEigenProblem
     and SolveEigenProblem, to be implemented in base classes.
     * @param[out] w Stores the eigenvalues
     * @return Returns 1 if executed correctly
     */
    inline int Solve(TPZVec<CTVar> &w) override{
        if(this->IsGeneralised()) return SolveGeneralisedEigenProblem(w);
        else return SolveEigenProblem(w);
    }
    //!Whether the solver is set for a generalised eigenvalue problem
    inline bool IsGeneralised() const{
        return fIsGeneralised;
    }
    //!Configure the solver to solve a generalised eigenvalue problem
    virtual void SetAsGeneralised(bool isGeneralised){
        fIsGeneralised = isGeneralised;
    }
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
    //! Resets Matrices
    void ResetMatrix() override
    {
        TPZAutoPointer<TPZMatrix<TVar>> newA, newB;
        fMatrixA = newA;
        fMatrixB = newB;
    }

    /** @name Pardiso*/
    /** @{*/
    //!Gets PARDISO control of Matrix A
    //!Returns nullptr if not applicable (non-sparse matrix or matrix not set)
    virtual TPZPardisoSolver<TVar> *GetPardisoControlA(){
        return GetPardisoControl(fMatrixA);
    }
    //!Gets PARDISO control of Matrix A
    //!Returns nullptr if not applicable (non-sparse matrix or matrix not set)
    virtual TPZPardisoSolver<TVar> *GetPardisoControlB(){
        return GetPardisoControl(fMatrixB);
    }
    /** @}*/
protected:
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
    /** @brief Whether to solve the eigenvalue problem
     *   is generalised (Ax=uBx) or not (Ax=ux)*/
    bool fIsGeneralised{false};
    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar>> fMatrixA{nullptr};

    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar>> fMatrixB{nullptr};

    TPZPardisoSolver<TVar> * GetPardisoControl(TPZAutoPointer<TPZMatrix<TVar>> mat);
};

extern template class TPZLinearEigenSolver<float>;
extern template class TPZLinearEigenSolver<double>;
extern template class TPZLinearEigenSolver<long double>;

extern template class TPZLinearEigenSolver<std::complex<float>>;
extern template class TPZLinearEigenSolver<std::complex<double>>;
extern template class TPZLinearEigenSolver<std::complex<long double>>;
#endif