/**
 * @file
 * @brief Contains TPZNonLinMultGridAnalysis class which implements multigrid analysis to non linear problems.
 */

#ifndef TPZNLMGANALYSIS_H
#define TPZNLMGANALYSIS_H

#include <iosfwd>                 // for string, ostream
#include "TPZLinearAnalysis.h"           // for TPZLinearAnalysis
#include "pzmatrix.h"             // for TPZFMatrix, TPZMatrix
#include "pzreal.h"               // for STATE, REAL
#include "pzstack.h"              // for TPZStack

class TPZCompMesh;
class TPZMaterial;
template <class TVar> class TPZMatrixSolver;


/**
 * @brief Implements multigrid analysis to Non linear problems. Class TPZNonLinMultGridAnalysis derived from TPZLinearAnalysis. \ref analysis "Analysis"
 * @ingroup Analysis
 */
class TPZNonLinMultGridAnalysis : public TPZLinearAnalysis {
	
	/** @brief Contains the computational meshes of one cycle: V, W, F, etc */
	TPZStack < TPZCompMesh * > fMeshes;
	
	/** @brief Contains the meshes solutions */	
	TPZStack <TPZFMatrix<STATE> *> fSolutions;
	
	/** @brief Contains the solution method applied to the mesh */
	TPZStack <TPZMatrixSolver<STATE> *> fSolvers;
	
	/** @brief Contains the preconditioner of the solution method */ 
	/** If the solution method is a krylov method, the preconditioner can be used as a coarse grid iteration */
	TPZStack <TPZMatrixSolver<STATE> *> fPrecondition;
	
	/** @brief Times by iteration and accumulated time */
	clock_t fBegin,fInit;
	
public:
	
	/** @brief Destructor */
	~TPZNonLinMultGridAnalysis();
	
	/** @brief Creates an object multigrid analysis giving a computational mesh */
	TPZNonLinMultGridAnalysis(TPZCompMesh *cmesh);
	
	/** @brief Append a mesh to the meshes vector */
	void AppendMesh (TPZCompMesh * mesh);
	
	/** @brief Pop the last mesh of the meshes vector */
	TPZCompMesh *PopMesh ();
	
	/** @brief Number of meshes */
	int NMeshes() {return fMeshes.NElements();}
	
	TPZCompMesh *IMesh(int64_t index);
	
	/**
	 * @brief It creates a new established computational mesh in the refinement uniform 
	 * of a fixed number of levels of the original geometric mesh
	 * @param coarcmesh : initial computational mesh
	 * @param levelnumbertorefine: number of the levels to be refined
	 * @param setdegree: degree of interpolation
	 * @note newmesh = 0 coarcmesh is refined other case new mesh is create
	 */
	static TPZCompMesh *UniformlyRefineMesh (TPZCompMesh *coarcmesh,int levelnumbertorefine,int setdegree);
	
	/** @brief It generates a new mesh based on the agglomeration of elements of the fine mesh */
	static TPZCompMesh *AgglomerateMesh(TPZCompMesh *finemesh,int levelnumbertogroup);
	
	void SmoothingSolution(REAL tol,int numiter,TPZMaterial * mat,TPZLinearAnalysis &an,int marcha = 0 ,
						   const std::string &dxout = "plotfile.dx");
	
	void SmoothingSolution(REAL tol,int numiter,TPZMaterial * mat,TPZLinearAnalysis &an,TPZFMatrix<STATE> &rhs);
	
	void SmoothingSolution2(REAL tol,int numiter,TPZMaterial * mat,TPZLinearAnalysis &an,int marcha,
							const std::string &dxout);
	
	void ResetReference(TPZCompMesh *aggcmesh);
	
	void SetReference(TPZCompMesh *aggcmesh);
	
	void SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial * mat);
	
	void CoutTime(clock_t &start,const char *title);
	
	void OneGridAlgorithm(std::ostream &out,int nummat);
	
	void TwoGridAlgorithm(std::ostream &out,int nummat);

	template<class TVar>
	void CalcResidual(TPZMatrix<TVar> &sol,TPZLinearAnalysis &an,const std::string  &decompose,TPZFMatrix<TVar> &res);

	template<class TVar>
	void CalcResidual(TPZMatrix<TVar> &sol,TPZFMatrix<TVar> &anres,TPZFMatrix<TVar> &res,TPZLinearAnalysis &an,const std::string &decompose);
		
public:
	
	void SetAnalysisFunction(void (*fp)(TPZMaterial *mat,TPZCompMesh *cmesh)){
		fFunction = fp;
	}
protected:
	
	void (*fFunction)(TPZMaterial *mat,TPZCompMesh *cmesh);
	
};

extern template void
TPZNonLinMultGridAnalysis::CalcResidual<STATE>(
	TPZMatrix<STATE> &sol, TPZLinearAnalysis &an,
	const std::string &decompose,
	TPZFMatrix<STATE> &res);

extern template void
TPZNonLinMultGridAnalysis::CalcResidual<STATE>(
    TPZMatrix<STATE> &sol, TPZFMatrix<STATE> &anres,
	TPZFMatrix<STATE> &res, TPZLinearAnalysis &an,
	const std::string &decompose);

#endif
