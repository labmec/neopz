/**
 * @file
 * @brief Contains TPZMGAnalysis class which implements multigrid analysis.
 */

#ifndef TPZMGANALYSIS_H
#define TPZMGANALYSIS_H

#include "pzanalysis.h"

class TPZTransfer;
class TPZInterpolatedElement;
class TPZTransform;
class TPZStepSolver;

template<class T, class V>

class TPZAvlMap;
class TPZOneDRef;
class TPZGeoEl;


/**
 * @brief Implements multigrid analysis. TPZMGAnalysis is derived from TPZAnalysis. \ref analysis "Analysis"
 * @ingroup Analysis
 */
class TPZMGAnalysis : public TPZAnalysis {
public:
	
	/**
	 * @brief Destructor
	 */
	virtual ~TPZMGAnalysis();
	
	/**
	 * @brief Creates an object multigrid analysis
	 * giving a computational mesh
	 */
	TPZMGAnalysis (TPZCompMesh *);
	
	/**
	 * @brief Append a mesh to the meshes vector
	 */
	void AppendMesh (TPZCompMesh * mesh);
	
	/**
	 * @brief Pop the last mesh of the meshes vector
	 */
	TPZCompMesh *PopMesh ();
	
	/**
	 * @brief Uses fSolver object to apply a solution
	 * algorithm
	 */
	virtual void Solve ();
	
	/**
	 * @brief Loads the last two solutions and
	 * call the error between these two aproximations
	 */
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
							void (*f) (TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),
							TPZVec<REAL> &truervec);
	
private:    
	/** @brief Contains the computational meshes of one cycle */
	TPZStack < TPZCompMesh * > fMeshes;
	
	/** @brief Contains the meshes solutions */	
	TPZStack <TPZFMatrix *> fSolutions;
	
	/** @brief Contains the solution method applied to the mesh */
	TPZStack <TPZMatrixSolver *> fSolvers;
	
	/** @brief Contains the preconditioner of the solution method if the solution method is a krylov method. */
	/** The preconditioner can be used as a coarse grid iteration */
	TPZStack <TPZMatrixSolver *> fPrecondition;
	
	
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
							   TPZTransform &tr,
							   void (*f) (TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),
							   REAL &truerror);
	
};

#endif //TPZMGANALYSIS_H
