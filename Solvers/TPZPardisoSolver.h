//  TPZPardisoSolver.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//
#ifndef TPZPARDISOSOLVER_H
#define TPZPARDISOSOLVER_H

#include "TPZMatrixSolver.h"

#include "pzmanvector.h"
#include "tpzautopointer.h"

template<class TVar>
class TPZSYsmpMatrixPardiso;

template<class TVar>
class TPZFYsmpMatrixPardiso;

/**
   @brief This class acts as a wrapper for calls to the PARDISO solver.
   It is both used internally by the sparse matrix classes or it can be used 
   by the user instead of TPZStepSolver to have a finer control over PARDISO parameters.
   Inspired by PardisoSolver by Armando Duarte.
   @ingroup solver
   @note It is only compatible with sparse matrices. Use it with TPZFYsmpMatrix or TPZSYsmpMatrix
   @note This class is meaningless unless the library is linked against MKL library.
*/
template<class TVar>
class TPZPardisoSolver : public TPZMatrixSolver<TVar>
{
    //they need access to SetRawMatrix and ReallySolve
    friend class TPZFYsmpMatrixPardiso<TVar>;
    friend class TPZSYsmpMatrixPardiso<TVar>;
public:
    /*! Whether the algebraic system is positive definite*/
    enum class MProperty {ENonInitialized=0, EPositiveDefinite, EIndefinite};
    
    /// Default constructor.
    TPZPardisoSolver();
    /*!Copy constructor*/
    TPZPardisoSolver(const TPZPardisoSolver &copy) = default;
    /*!Move constructor*/
    TPZPardisoSolver(TPZPardisoSolver &&copy) = default;
    /*!Copy-assignment operator*/
    TPZPardisoSolver &operator=(const TPZPardisoSolver &copy) = default;
    /*!Move-assignment operator*/
    TPZPardisoSolver &operator=(TPZPardisoSolver &&copy) = default;
    /*!Destructor*/
    virtual ~TPZPardisoSolver();
    /** @brief Sets the associated matrix of the solver.
       If SetMatrixType has not been called, the following defaults are assumed:
       Symmetric matrices: `MSystemType::ESymmetric`, `MProperty::EIndefinite` 
       Non-symmetric matrices: `MSystemType::ENonSymmetric`, `MProperty::EIndefinite`
       @note Only sparse matrices are compatible with this solver.*/
    void SetMatrix(TPZAutoPointer<TPZBaseMatrix> Refmat) override;
    /**@brief Decompose the matrix */
    void Decompose() override;
    
    /** @brief Use the decomposed matrix to invert the system of equations
     This method also decomposed the matrix if needed.*/
    void Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol, TPZFMatrix<TVar> *residual=nullptr) override;

    /**
       @brief Change the matrix type
       @note  This method should only be called by the user before a matrix has been set.
    */
    void SetMatrixType(SymProp symtype, MProperty prop);

    //! Changes only PARDISO's matrix type param
    void SetMatrixType(long long type) {fMatrixType = type;} 
    /** @brief Sets the message level of the pardiso solver.
     Anything different than zero results in statistical information being printed.*/
    void SetMessageLevel(int lvl);
    //! Clones the current object returning a pointer of type TPZSolver
	  TPZPardisoSolver<TVar> *Clone() const override;
    /** @brief Gets copy of Pardiso param array
        @note It is strongly suggested to call SetStructure and SetMatrixType before
        customising the param array.*/
    [[nodiscard]] inline TPZVec<long long> GetParam() const
    {return fParam;}

    /** @brief Sets custom values to Pardiso param array.*/
    void SetParam(const TPZVec<long long> & p);

    /** @brief Resets custom values of Pardiso param array
     @note This function should only be called if SetParam has been previously called.*/
    void ResetParam();
    //! Whether the param array has been set manually.
    [[nodiscard]] bool HasCustomSettings() const {return fCustomSettings;}

    //! Returns permutation vec used by pardiso (can be obtained by param[4]=2)
    TPZVec<long long> &GetPermutationVec() {return fPermutation;}

    //! Releases all internal memory from PARDISO solver
    void FreePardisoMemory();
protected:
    /// Compute the `mtype` parameter of the pardiso_64 call
    long long MatrixType();
    /**@brief Decompose the matrix */
    void Decompose(TPZMatrix<TVar> *mat);
    /** @brief Use the decomposed matrix to invert the system of equations
     This method assumes that the matrix has been decomposed.*/
    void Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const;
    
    SymProp fSymmetry{SymProp::NonSym};
    
    MProperty fProperty{MProperty::ENonInitialized};
    
    // Solver internal data address pointers
    // 32-bit: int pt[64]; 64-bit: long long pt[64]
    // or void *pt[64] should be OK on both architectures
    // this datastructure should not be copied or duplicated, therefore the "autopointer" protection
    TPZAutoPointer<TPZManVector<long long, 64> > fPardisoControl{nullptr};
    
    // adress of the first element of pt;
    long long *fHandle{nullptr};
    
    // Array used to pass parameters to Pardiso
    TPZManVector<long long, 64> fParam;
    
    // Maximum number of factors we will pass to the solver
    long long fMax_num_factors{1};
    
    // Factor number we are using
    long long fMatrix_num{1};
    
    // Message level information
    long long fMessageLevel{0};
    
    // error flag from Pardiso
    long long fError{0};
    
    /// permutation vector computed by Pardiso
    TPZVec<long long> fPermutation;
    
    // matrix type, computed based on the structural information and TVar
    long long fMatrixType{0};

    /// whether the matrix has been decomposed
    bool fDecomposed{false};
    /// whether pardisoinit has been called
    bool fPardisoInitialized{false};
    /// whether a custom param array has been set
    bool fCustomSettings{false};
};

#endif /* TPZPARDISOSOLVER_H */
