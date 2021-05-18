/**
 * @file
 * @brief Contains TPZMatrixSolver class which defines an abstract class of solvers which will be used by matrix classes for solving equations systems.
 */

#ifndef TPZMATRIXSOLVER_H
#define TPZMATRIXSOLVER_H

#include "pzfmatrix.h"
#include "TPZSolver.h"
#include "tpzautopointer.h"



/**
 * @ingroup solver
 * @brief  Defines the interface for matrix solvers. \ref solver "Solver"
 * These are used for solving linear or non-linear systems of equations.
 */
template<class TVar>
class TPZMatrixSolver: public TPZSolver
{
	
public:
    /*! Defines the solvers available in NeoPZ*/
	enum MSolver
	{
		ENoSolver, /*!< No solver selected*/
        EJacobi,/*!< Jacobi solver*/
        ESOR,/*!<Successive Over Relaxation solver*/
        ESSOR,/*!< Symmetric Successive Over Relaxation solver*/
        ECG,/*!< Conjugate Gradient solver*/
        EGMRES,/*!<Generalised Minimal Residual method*/
        EBICGSTAB,/*!<Biconjugate Gradient Stabilized Method*/
        EDirect,/*!< LU, LDLt or Cholesky*/
        EMultiply
	};

	/**
	 * @brief Constructor with initialization parameter
	 * @param Refmat Sets reference matrix to `nullptr`
	 */
	
	TPZMatrixSolver(TPZAutoPointer<TPZMatrix<TVar> >  Refmat);
	
	TPZMatrixSolver();
	
	/**
	 * @brief Copy constructor
	 * @param Source Model object to be copied from
	 */
	TPZMatrixSolver(const TPZMatrixSolver<TVar> &Source);
	
	/** @brief Destructor */
	virtual ~TPZMatrixSolver();

    /**
     * @brief Solves the system of linear equations
     * @param F contains Force vector
     * @param result contains the solution
     * @param residual contains the residual for that linear system
     */
    virtual void Solve(const TPZFMatrix<TVar> &F,
                       TPZFMatrix<TVar> &result,
                       TPZFMatrix<TVar> *residual = 0) = 0;
    
    /** @brief Decompose the system of equations if a direct solver is used*/
    virtual void Decompose() {}

        /**
	 * @brief Sets a matrix to the current object
	 * @param Refmat Sets reference matrix to RefMat
	 */
	virtual void SetMatrix(TPZAutoPointer<TPZBaseMatrix> Refmat)
    {
        fContainer = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(Refmat);
    }
	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZBaseMatrix> mat_base)
	{
        auto matrix =
            TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(mat_base);
		if (fReferenceMatrix == matrix && matrix)
		{
			if(this->fContainer) this->fContainer->UpdateFrom(matrix);
		}
	}
	/** @brief Resets current object */
	void ResetMatrix() override;
	
	/** @brief This method gives a preconditioner to share a matrix with the referring solver object */
	virtual void SetReferenceMatrix(TPZAutoPointer<TPZBaseMatrix> matrix)
	{
		fReferenceMatrix = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(matrix);
	}
	
	/** @brief Returns a pointer to TPZMatrix<>*/
	TPZAutoPointer<TPZMatrix<TVar> > Matrix() const
	{
		return fContainer;
	}
	
    void ReallocMatrix() {
        fContainer.ReallocForNuma(0);
    }
    
	/**
	 * @brief Shares the current matrix with another object of same type
	 * @param other Object that will share current matrix
	 */
	void ShareMatrix(TPZMatrixSolver<TVar> & other);
	
    virtual MSolver Solver()
    {
        return ENoSolver;
    }
    /** @brief Saveable specific methods */
    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;

protected:
    /** @brief Reference matrix used to update the current matrix */
    TPZAutoPointer<TPZMatrix<TVar>> fReferenceMatrix;

    /** @brief Manipulation matrix */
    TPZFMatrix<TVar> fScratch;

private:
    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar>> fContainer;	
};

extern template class TPZMatrixSolver<float>;
extern template class TPZMatrixSolver<std::complex<float> >;

extern template class TPZMatrixSolver<double>;
extern template class TPZMatrixSolver<std::complex<double> >;

extern template class TPZMatrixSolver<long double>;
extern template class TPZMatrixSolver<std::complex<long double> >;

extern template class TPZMatrixSolver<Fad<float> >;
extern template class TPZMatrixSolver<Fad<double> >;
extern template class TPZMatrixSolver<Fad<long double> >;
#endif  // TPZMATRIXSOLVER_H
