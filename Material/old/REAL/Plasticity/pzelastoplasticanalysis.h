/**
 * @file
 * @brief Contains the declaration of TPZElastoPlasticAnalysis class
 */

#ifndef ELASTOPLASTICANALYSIS_H
#define ELASTOPLASTICANALYSIS_H

#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "TPZGeoElement.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzpostprocanalysis.h"
#include <iostream>


class TPZElastoPlasticAnalysis : public TPZNonLinearAnalysis {

public:
	/** @brief Constructor */
	TPZElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out);
	/** @brief Default constructor */
	TPZElastoPlasticAnalysis();
	/** @brief Default destructor */
	virtual ~TPZElastoPlasticAnalysis();
	
    virtual void IterativeProcessPrecomputedMatrix(std::ostream &out, REAL tol, int numiter, bool linesearch);
    
    /// Iterative process using the linear elastic material as tangent matrix
    virtual void IterativeProcess(std::ostream &out, REAL tol, int numiter, int niter_update_jac, bool linesearch);
    
	/**
	 * @brief It process a Newton's method to solve the non-linear problem.
	 * In this implementation, the line search is temporarily disabled.
	 */
	virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,bool &ConvOrDiverg);

    	virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv);
    
	//virtual REAL LineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter);
    
    //Improved:A verification is made in order to check convergence.
    virtual REAL LineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool &converging );
	
	/**
	 * @brief The code below manages the update of a certain boundary condition (BCId)
	 * to assume values progressing from begin to end values within nsteps, such
	 * that the last step length is a 'lastStepRatio' fraction of the global
	 * load in a geometric progression.
	 * @param out [in] output device to write the output to
	 * @param tol [in] tolerance desired in each loading step
	 * @param numiter [in] maximum iteration steps in each loading step
	 * @param BCId [in] boundary condition Id
	 * @param nsteps [in] number of steps
	 * @param PGRatio [in] ratio for the PG progression
	 * @param val1Begin [in] beginning value of BC.Val1
	 * @param val1End [in] final value of BC.Val1
	 * @param val2Begin [in] beginning value of BC.Val2
	 * @param val2End [in] final value of BC.Val2
	 * @param ppAnalysis [in] TPZPostProcAnalysis object to write the output to
	 * @param res [in] output postprocessing resolution
	 */
	virtual void ManageIterativeProcess(std::ostream &out,REAL tol,int numiter,
										int BCId, int nsteps, REAL PGRatio,
										TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
										TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
										TPZPostProcAnalysis * ppAnalysis = NULL, int res = 0);
	
	/**
	 * @brief Informs the materials to update the plastic memory, assembles the
	 * rhs in order to update the memory and unsets the materials
	 * to update the memory.
	 * @param ResetOutputDisplacements [in] Informs whether to add or reset the deltaSolution to the cumulative solution.
	 */
	virtual REAL AcceptSolution(const int ResetOutputDisplacements = 0);
    
    /** @brief Load the solution into the computable grid, transfering it to the multi physics meshes */
    virtual void LoadSolution();

    
    TPZFMatrix<REAL> &CumulativeSolution()
    {
        return fCumSol;
    }
	
	void SetPrecond(TPZMatrixSolver<REAL> &precond);
	
	void SetBiCGStab(int numiter, REAL tol);
	
	void SetBiCGStab_Jacobi(int numiter, REAL tol);
	
	void SetLU();
	
	void TransferSolution(TPZPostProcAnalysis & ppanalysis);
	
    void AddNoPenetration(int matid, int direction)
    {
        fMaterialIds.insert(std::pair<int,int>(matid, direction));
    }
    
    /// build the fEquationstoZero datastructure based on the fMaterialIds data structure
    void IdentifyEquationsToZero();
    
    /// return the vector of active equation indices
    void GetActiveEquations(TPZVec<int64_t> &activeEquations);
    
protected:	
	
	/**
	 * @brief Forces the materials with memory to update the internal
	 * plastic memory during the subsequent assemble/contribute calls
	 * when update is set to 1. When set to 0, the conventional
	 * contribute routine takes place
	 * @param update [in] 1 - updates the material memory during the next assemble calls; 0 - conventional assemble
	 */
	void SetUpdateMem(int update);
	
	/**
	 * @brief Updates block diagonal preconditioning matrix.
	 */
	void UpdatePrecond();
	
	
public:
	
	void CheckConv(std::ostream &out, REAL range);
	
	virtual void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase);
	
	virtual int NumCases();
	
	virtual void Residual(TPZFMatrix<REAL> &residual, int icase);
	
	static void SetAllCreateFunctionsWithMem(TPZCompMesh *cmesh);
	
    /// Structure defining a multiphysics configuration
    bool IsMultiPhysicsConfiguration()
    {
        return fMultiPhysics != NULL;
    }
    
    /// Initialize the multiphysics data structure
    void SetMultiPhysics(TPZCompMesh *mphysics, TPZVec<TPZCompMesh *> &meshvec)
    {
        fMultiPhysics = mphysics;
        fMeshVec = meshvec;
    }
    
    /// Reset the multiphysics data structure
    void ResetMultiPhysics()
    {
        fMultiPhysics = 0;
        fMeshVec.Resize(0);
    }
    
	
protected:
	
	/* @brief Cumulative solution vector*/
	TPZFMatrix<REAL> fCumSol;
	
	TPZMatrixSolver<REAL> * fPrecond;
    
    /// Equations with zero dirichlet boundary condition
    std::set<int64_t> fEquationstoZero;
    
    /// Materials with no penetration boundary conditions
    // the second value of the map indicates x (0) or y (1) restraint
    std::multimap<int,int> fMaterialIds;
    
    /// The multiphysics mesh
    TPZCompMesh *fMultiPhysics;
    
    /// The elasticity mesh and vertical deformation mesh
    TPZManVector<TPZCompMesh *,2> fMeshVec;
	
	/**
	 * TPZCompElWithMem<TBASE> creation function setup
	 * These functions should be defined as static members of the TPZCompElWithMem<TBASE> class
	 * but were defined here because the TPZCompElWithMem is a template class and would
	 * require a dummy template argumet in order to be called. That woudn't be elegant.
	 */
	
	static TPZCompEl * CreateCubeElWithMem(  TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateLinearElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreatePointElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreatePrismElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreatePyramElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateQuadElWithMem(  TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateTetraElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateTriangElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	
};

#endif

