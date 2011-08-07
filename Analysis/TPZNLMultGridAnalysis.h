/**
 * @file
 * @brief Contains TPZNonLinMultGridAnalysis class which implements multigrid analysis to non linear problems.
 */

#ifndef TPZNLMGANALYSIS_H
#define TPZNLMGANALYSIS_H

#include "pzanalysis.h"

class TPZTransfe;
class TPZInterpolatedElement;
class TPZTransform;
class TPZStepSolver;
class TPZCompMesh;

template<class T, class V>

class TPZAvlMap;
class TPZOneDRef;
class TPZGeoEl;


/**
 * @brief Implements multigrid analysis to Non linear problems. Class TPZNonLinMultGridAnalysis derived from TPZAnalysis. \ref analysis "Analysis"
 * @ingroup Analysis
 */
class TPZNonLinMultGridAnalysis : public TPZAnalysis {
	
	/**
	 * @brief Contains the computational meshes of one cycle: V, W, F, etc
	 */
	TPZStack < TPZCompMesh * > fMeshes;
	
	/**
	 * @brief Contains the meshes solutions
	 */	
	TPZStack <TPZFMatrix *> fSolutions;
	
	/**
	 * @brief Contains the solution method applied to the mesh
	 */
	TPZStack <TPZMatrixSolver *> fSolvers;
	
	/**
	 * @brief Contains the preconditioner of the solution method
	 * 
	 * If the solution method is a krylov method, the preconditioner
	 * can be used as a coarse grid iteration
	 */
	TPZStack <TPZMatrixSolver *> fPrecondition;
	
	clock_t fBegin,fInit;
	
public:
	
	/**
	 * @brief Destructor
	 */
	~TPZNonLinMultGridAnalysis();
	
	/**
	 * @brief Creates an object multigrid analysis
	 * giving a computational mesh
	 */
	TPZNonLinMultGridAnalysis(TPZCompMesh *cmesh);
	
	/**
	 * @brief Append a mesh to the meshes vector
	 */
	void AppendMesh (TPZCompMesh * mesh);
	
	/**
	 * @brief Pop the last mesh of the meshes vector
	 */
	TPZCompMesh *PopMesh ();
	
	int NMeshes() {return fMeshes.NElements();}
	
	TPZCompMesh *IMesh(int index);
	
	/**
	 * Uses fSolver object to apply a solution
	 * algorithm
	 */
	//virtual void Solve ();
	
	/**
	 * @brief It creates a new established computational mesh in the refinement uniform 
	 * of a fixed number of levels of the original geometric mesh
	 * @param coarcmesh : initial computational mesh
	 * @param levelnumbertorefine: number of the levels to be refined
	 * @param setdegree: degree of interpolation
	 * @note newmesh = 0 coarcmesh is refined other case new mesh is create
	 */
	static TPZCompMesh *UniformlyRefineMesh (TPZCompMesh *coarcmesh,int levelnumbertorefine,int setdegree);
	
	/**
	 * @brief It generates a new mesh based on the agglomeration of elements of the fine mesh
	 */
	static TPZCompMesh *AgglomerateMesh(TPZCompMesh *finemesh,int levelnumbertogroup);
	
	void SmoothingSolution(REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,TPZAnalysis &an,int marcha = 0 ,
						   const std::string &dxout = "plotfile.dx");
	
	void SmoothingSolution(REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,TPZAnalysis &an,TPZFMatrix &rhs);
	
	void SmoothingSolution2(REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,TPZAnalysis &an,int marcha,
							const std::string &dxout);
	
	void ResetReference(TPZCompMesh *aggcmesh);
	
	void SetReference(TPZCompMesh *aggcmesh);
	
	void SetDeltaTime(TPZCompMesh *CompMesh,TPZAutoPointer<TPZMaterial> mat);
	
	void CoutTime(clock_t &start,const char *title);
	
	void OneGridAlgorithm(std::ostream &out,int nummat);
	
	void TwoGridAlgorithm(std::ostream &out,int nummat);
	
	void CalcResidual(TPZMatrix &sol,TPZAnalysis &an,const std::string  &decompose,TPZFMatrix &res);
	
	void CalcResidual(TPZMatrix &sol,TPZFMatrix &anres,TPZFMatrix &res,TPZAnalysis &an,const std::string &decompose);
	
	/*   void IterativeProcess(TPZAnalysis &an,REAL tol,int numiter, */
	/* 			TPZMaterial *mat,int marcha,int resolution); */
	
	
protected:
	
	void (*fFunction)(TPZMaterial *mat,TPZCompMesh *cmesh);
	
public:
	
	void SetAnalysisFunction(void (*fp)(TPZMaterial *mat,TPZCompMesh *cmesh)){
		fFunction = fp;
	}
	
};

#endif
