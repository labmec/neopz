/**
 * @file
 * @brief Contains TPZSolver class which defines a abstract class of solvers  which will be used by matrix classes.
 */

#ifndef TPREH
#define TPREH

#include "pzfmatrix.h"
#include "TPZHash.h"
template<class TVar>
class TPZMatrixSolver;

/**
 * @ingroup solver
 * @brief Defines a abstract class of solvers  which will be used by matrix classes. \ref solver "Solver"
 */
template<class TVar>
class TPZSolver: public TPZSavable
{

public:
    
    public:
int ClassId() const override;

	/**
	 * @brief Solves the system of linear equations
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual contains the residual for that linear system
	 */
	virtual void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
					   TPZFMatrix<TVar>  *residual = 0) = 0;
    
    /** @brief Decompose the system of equations if a direct solver is used */
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
	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		std::cout << __PRETTY_FUNCTION__ << " called\n";
	}

};


template<class TVar>
int TPZSolver<TVar>::ClassId() const{
    return Hash("TPZSolver") ^ ClassIdOrHash<TVar>() << 1;
}

/** @ingroup solver */
#define TPZMATRIXSOLVER_ID 28291005;


/**
 * @ingroup solver
 * @brief  Defines a class of matrix solvers. \ref solver "Solver"
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
	 * @brief Constructor with initialization parameter
	 * @param Refmat Sets reference matrix to 0
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
	 * @brief Sets a matrix to the current object
	 * @param Refmat Sets reference matrix to RefMat
	 */
	virtual void SetMatrix(TPZAutoPointer<TPZMatrix<TVar> > Refmat)
    {
        fContainer = Refmat;
    }
	
	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix) override
	{
		if (fReferenceMatrix == matrix && matrix)
		{
			if(this->fContainer) this->fContainer->UpdateFrom(matrix);
		}
	}
	/** @brief Resets current object */
	void ResetMatrix() override;
	
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

protected:
	
private:
	/** @brief Container classes */
	TPZAutoPointer<TPZMatrix<TVar> > fContainer;
protected:
	/** @brief Reference matrix used to update the current matrix */
	TPZAutoPointer<TPZMatrix<TVar> > fReferenceMatrix;

protected:
	/** @brief Manipulation matrix */
	TPZFMatrix<TVar>  fScratch;
public:
	/** @brief Saveable specific methods */
	public:
int ClassId() const override;

	void Write(TPZStream &buf, int withclassid) const override;
	void Read(TPZStream &buf, void *context) override;
};

template<class TVar>
int TPZMatrixSolver<TVar>::ClassId() const{
    return Hash("TPZMatrixSolver") ^ TPZSolver<TVar>::ClassId() << 1;
}

#endif  // TPREH
