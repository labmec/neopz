/**
 * @file
 * @brief Contains TPZSolver class which defines a abstract class of solvers  which will be used by matrix classes.
 */
#ifndef TPREH
#define TPREH

#include "pzfmatrix.h"

template<class TVar>
class TPZMatrixSolver;

/**
 @ingroup solver
 @brief Defines a abstract class of solvers  which will be used by matrix classes. \ref solver "Solver"
 */
template<class TVar>
class TPZSolver: public TPZSaveable
{

public:
	/**
	 * @brief Solves the system of linear equations
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual contains the residual for that linear system
	 */
	virtual void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
					   TPZFMatrix<TVar>  *residual = 0) = 0;
    
    /**
     * @brief Decompose the system of equations if a direct solver is used
     */
    virtual void Decompose()
    {
    }
	
	/** @brief Clones the current object returning a pointer of type TPZSolver */
	virtual TPZSolver *Clone() const = 0;
	/** @brief Destructor */
	virtual ~TPZSolver();
	
	/** @brief This method will reset the matrix associated with the solver */
	/** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix()
	{
	}
	/**
	 This method gives a preconditioner to share a matrix with the referring solver object
	 */
	//  virtual void SetMatrix(TPZMatrixSolver *solver);
	/**
	 * @brief Updates the values of the current matrix based on the values of the matrix
	 */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		std::cout << __PRETTY_FUNCTION__ << " called\n";
	}

};

/** @ingroup solver */
#define TPZMATRIXSOLVER_ID 28291005;


/**
 @ingroup solver
 @brief  Defines a class of matrix solvers. \ref solver "Solver"
 */
template<class TVar>
class TPZMatrixSolver: public TPZSolver<TVar>
{
	
public:
    /**
	 * @enum MSolver
	 * @brief Defines a series of solvers available in PZ
	 * @param ENoSolver No solver selected
	 * @param EJacobi Jacobi solver selected
	 * @param ESOR Successive Over Relaxation solver selected
	 * @param ESSOR Symmetric Successive Over Relaxation solver selected
	 * @param ECG Conjugate Gradiente solver selected
	 * @param EDirect LU, LDLt or Cholesky selected
	 */
	enum MSolver
	{
		ENoSolver, EJacobi, ESOR, ESSOR, ECG, EGMRES, EBICGSTAB, EDirect, EMultiply
	};

	/**
	 @brief Constructor with initialization parameter
	 @param Refmat Sets reference matrix to 0
	 */
	
	TPZMatrixSolver(TPZAutoPointer<TPZMatrix<TVar> >  Refmat);
	
	TPZMatrixSolver();
	
	/**
	 @brief Copy constructor
	 @param Source Model object to be copied from
	 */
	TPZMatrixSolver(const TPZMatrixSolver<TVar> &Source);
	
	/** @brief Destructor */
	virtual ~TPZMatrixSolver();
	
	/**
	 @brief Sets a matrix to the current object
	 @param Refmat Sets reference matrix to RefMat
	 */
	virtual void SetMatrix(TPZAutoPointer<TPZMatrix<TVar> > Refmat);
	
	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		if (fReferenceMatrix == matrix && matrix)
		{
			if(this->fContainer) this->fContainer->UpdateFrom(matrix);
		}
	}
	/** @brief Resets current object */
	void ResetMatrix();
	
	/** @brief This method gives a preconditioner to share a matrix with the referring solver object */
	virtual void SetReferenceMatrix(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		fReferenceMatrix = matrix;
	}
	
	/** @brief Returns a pointer to TPZMatrix<>*/
	TPZAutoPointer<TPZMatrix<TVar> > Matrix() const
	{
		return fContainer;
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

protected:
	
private:
	/** @brief Container classes */
	TPZAutoPointer<TPZMatrix<TVar> > fContainer;
protected:
	/** @brief Reference matrix used to update the current matrix */
	TPZAutoPointer<TPZMatrix<TVar> > fReferenceMatrix;
	//	TPZSolver *fPrecond;

protected:
	/** @brief Manipulation matrix */
	TPZFMatrix<TVar>  fScratch;
public:
	/** @brief Saveable specific methods */
	virtual int ClassId() const
	{
		return TPZMATRIXSOLVER_ID;
	}
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
};

#endif  // TPREH
