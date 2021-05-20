/**
 * @file
 * @brief Contains TPZMGAnalysis class which implements multigrid analysis.
 */

#ifndef TPZMGANALYSIS_H
#define TPZMGANALYSIS_H

#include "TPZLinearAnalysis.h"
#include "pztransfer.h"
#include "pztrnsform.h" //needed because of default templ params

class TPZInterpolatedElement;
template<class T>
class TPZTransform;

template <class TVar>
class TPZStepSolver;

template<class T, class V>

class TPZAvlMap;
class TPZOneDRef;
class TPZGeoEl;

/**
 * @brief Implements multigrid analysis. TPZMGAnalysis is derived from TPZLinearAnalysis. \ref analysis "Analysis"
 * @ingroup Analysis
 */
class TPZMGAnalysis : public TPZLinearAnalysis {
public:
	
	/** @brief Destructor */
	virtual ~TPZMGAnalysis();
	
	/** @brief Creates an object multigrid analysis giving a computational mesh */
	TPZMGAnalysis (TPZCompMesh *);
	
	/** @brief Append a mesh to the meshes vector */
	void AppendMesh (TPZCompMesh * mesh);
	
	/** @brief Pop the last mesh of the meshes vector */
	TPZCompMesh *PopMesh ();
	
	/** @brief Uses fSolver object to apply a solution algorithm */
	virtual void Solve () override;
	
	/** @brief Loads the last two solutions and call the error between these two aproximations */
	void ComputeError (TPZVec<REAL> &error);
	
	/**
	 * @brief Proceeds the uniformly h-p refinement of mesh
	 * @param mesh : input mesh which will be refined
	 * @param withP : if true, increase the p-order
	 */
	static TPZCompMesh *UniformlyRefineMesh (TPZCompMesh *mesh, bool withP = false);
	
	/**
	 * @brief Evaluates the error between aproximation of coarse and fine meshes
	 * @param fine refined mesh
	 * @param coarse some father mesh of fine
	 * @param ervec	will return the calculated element error
	 * @param f pointer function
	 * @param truervec calculates the true error between a giving a solution
	 */
	static  void MeshError (  TPZCompMesh *fine,TPZCompMesh *coarse,	
							TPZVec<REAL> &ervec,
							void (*f) (const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),
							TPZVec<REAL> &truervec);
	
private:    
	/** @brief Contains the computational meshes of one cycle */
	TPZStack < TPZCompMesh * > fMeshes;
	
	/** @brief Contains the meshes solutions */	
	TPZStack <TPZSolutionMatrix *> fSolutions;
	
	/** @brief Contains the solution method applied to the mesh */
	TPZStack <TPZMatrixSolver<STATE> *> fSolvers;
	
	/** @brief Contains the preconditioner of the solution method if the solution method is a krylov method. */
	/** The preconditioner can be used as a coarse grid iteration */
	TPZStack <TPZMatrixSolver<STATE> *> fPrecondition;
	
	
	/**
	 * @brief Calculates an element error based on two aproximations
	 * @param fine refined mesh;
	 * @param coarse some father mesh of fine;
	 * @param tr tranformation between fine and coarse;
	 * @param f pointer function
	 * @param truerror will return the error between aproximation and a give solution
	 */
	static  REAL ElementError (TPZInterpolatedElement *fine,
							   TPZInterpolatedElement *coarse,
							   TPZTransform<> &tr,
							   void (*f) (const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),
							   REAL &truerror);
	template<class TVar>
	void SolveInternal();
};

#endif //TPZMGANALYSIS_H
