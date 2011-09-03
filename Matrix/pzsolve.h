/**
 * @file
 * @brief Contains TPZSolver class which defines a abstract class of solvers  which will be used by matrix classes.
 */
#ifndef TPREH
#define TPREH

#include "pzfmatrix.h"

class TPZMatrixSolver;

/**
 @ingroup solver
 @brief Defines a abstract class of solvers  which will be used by matrix classes. \ref solver "Solver"
 */
class TPZSolver: public TPZSaveable
{

public:
	// virtual void Solve(TPZFMatrix &F, TPZFMatrix &result) = 0;
	/**
	 * @brief Solves the system of linear equations stored in current matrix
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual contains the residual for that linear system
	 */
	virtual void Solve(const TPZFMatrix &F, TPZFMatrix &result,
					   TPZFMatrix *residual = 0) = 0;
	
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
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix> matrix)
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
class TPZMatrixSolver: public TPZSolver
{
	
public:
	/**
	 @brief Constructor with initialization parameter
	 @param Refmat Sets reference matrix to 0
	 */
	TPZMatrixSolver(TPZAutoPointer<TPZMatrix> Refmat);
	
	TPZMatrixSolver();
	
	/**
	 @brief Copy constructor
	 @param Source Model object to be copied from
	 */
	TPZMatrixSolver(const TPZMatrixSolver &Source);
	
	/** @brief Destructor */
	virtual ~TPZMatrixSolver();
	
	/**
	 @brief Sets a matrix to the current object
	 @param Refmat Sets reference matrix to RefMat
	 */
	virtual void SetMatrix(TPZAutoPointer<TPZMatrix> Refmat);
	
	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix> matrix)
	{
		if (fReferenceMatrix == matrix && matrix)
		{
			if(this->fContainer) this->fContainer->UpdateFrom(matrix);
		}
	}
	/** @brief Resets current object */
	void ResetMatrix();
	
	/** @brief This method gives a preconditioner to share a matrix with the referring solver object */
	virtual void SetReferenceMatrix(TPZAutoPointer<TPZMatrix> matrix)
	{
		fReferenceMatrix = matrix;
	}
	
	/** @brief Returns a pointer to TPZMatrix */
	TPZAutoPointer<TPZMatrix> Matrix() const
	{
		return fContainer;
	}
	
	/**
	 * @brief Shares the current matrix with another object of same type
	 * @param other Object that will share current matrix
	 */
	void ShareMatrix(TPZMatrixSolver & other);
	
protected:
	/**
	 * @enum MSolver
	 * @brief Defines a series of solvers available in PZ
	 * @param ENoSolver No solver selected
	 * @param EJacobi Jacobi solver selected
	 * @param ESOR Successive Over Relaxation solver selected
	 * @param ESSOR Symmetric Successive Over Relaxation solver selected
	 * @param ECG Conjugate Gradiente solver selected
	 * @param EDirect Jacobi solver selected
	 */
	enum MSolver
	{
		ENoSolver, EJacobi, ESOR, ESSOR, ECG, EGMRES, EBICGSTAB, EDirect, EMultiply
	};
	
private:
	/**
	 * @brief Defines a class of containers. \ref solver "Solver"
	 * @ingroup solver
	 */
	class TPZContainer
	{
		static int gnumcreated;
		static int gnumdeleted;
		int fRefCount;
		TPZMatrix *fRefMat;
	public:
		TPZContainer(TPZMatrix *mat);
		~TPZContainer();
		void IncreaseRefCount();
		void DecreaseRefCount();
		TPZMatrix *Matrix();
		void SetMatrix(TPZMatrix *mat);
		static void Diagnose(std::ostream &out = std::cout);
	};
	
	/** @brief Container classes */
	TPZAutoPointer<TPZMatrix> fContainer;
protected:
	/** @brief Reference matrix used to update the current matrix */
	TPZAutoPointer<TPZMatrix> fReferenceMatrix;
	//	TPZSolver *fPrecond;

protected:
	/** @brief Manipulation matrix */
	TPZFMatrix fScratch;
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
