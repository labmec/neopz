/**
 * @file
 * @brief Contains TPZAnalysis class which implements the sequence of actions to perform a finite element analysis.
 */
#ifndef ANALYSISH
#define ANALYSISH

class TPZGeoMesh;
class TPZCompMesh;
template<class TVar>
class TPZMatrix;
template<class TVar>
class TPZBlock;
class TPZConnect;
template<class TVar>
class TPZSolver;
template<class TVar>
class TPZMatrixSolver;
class TPZCompEl;
class TPZGraphMesh;
class TPZMaterial;
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzrenumbering.h"
#include "pzstrmatrix.h"
#include <iostream>


#include "pzfmatrix.h"
#include "TPZGuiInterface.h"

class TPZStructMatrix;
template<class T, int N> class TPZStack;

/**
 * @ingroup analysis
 * @brief Implements the sequence of actions to perform a finite element analysis. \ref analysis "Analysis"
 */
/** This class will renumerate the nodes upon construction
 */
class TPZAnalysis {
	
public:
	
	/** @brief  Preconditioners which can be created by objects of this class */
	enum EPrecond { EJacobi, EBlockJacobi, EElement, ENodeCentered };
	
	
protected:
	/**
	 * @brief Geometric Mesh
	 */
	TPZGeoMesh *fGeoMesh;
	/**
	 * @brief Computational mesh
	 */
	TPZCompMesh *fCompMesh;
	/**
	 * @brief Graphical mesh
	 */
	TPZGraphMesh *fGraphMesh[3];
	/**
	 * @brief Load vector ???
	 */
	TPZFMatrix<REAL> fRhs;
	/**
	 * @brief Solution vector
	 */
	TPZFMatrix<REAL> fSolution;
	/**
	 * @brief Type of solver to be applied
	 */
	TPZMatrixSolver<REAL> *fSolver;
	/**
	 * @brief Scalar variables names - to post process
	 */
	TPZVec<std::string> fScalarNames[3];
	/**
	 * @brief Vector variables names - to post process
	 */
	TPZVec<std::string> fVectorNames[3];
	/**
	 * @brief Step ???
	 */
	int fStep;
	/**
	 * @brief Time variable which is used in dx output
	 */
	REAL fTime;
	
	/**
	 * @brief Structural matrix
	 */
	TPZAutoPointer<TPZStructMatrix>  fStructMatrix;
	
	/**
	 * @brief Renumbering scheme
	 */
	TPZAutoPointer<TPZRenumbering> fRenumber;
	
	/** @brief Pointer for gui interface object */
	TPZAutoPointer<TPZGuiInterface> fGuiInterface;
	
	/** @brief Datastructure which defines postprocessing for one dimensional meshes */
	struct TTablePostProcess {
		TPZVec<int> fGeoElId;
		TPZVec<TPZCompEl *> fCompElPtr;
		int fDimension;
		TPZVec<REAL> fLocations;
		TPZVec<char *> fVariableNames;
		std::ostream *fOutfile;
		TTablePostProcess();
		~TTablePostProcess();
	} fTable;
	
	public :
	
    /**
     * @brief Exact solution ??? not implemented???
     */
    void (*fExact)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix<REAL> &deriv);
	
	/**
	 * @brief Create an TPZAnalysis object from one mesh pointer
	 */
	TPZAnalysis(TPZCompMesh *mesh,std::ostream &out = std::cout);
	
	/**
	 * @brief Create an TPZAnalysis object from one mesh pointer
	 */
	TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh,std::ostream &out = std::cout);
	
	/** @brief Defines gui interface object */
	void SetGuiInterface(TPZAutoPointer<TPZGuiInterface> gui){
		fGuiInterface = gui;
	}
	
	TPZAutoPointer<TPZGuiInterface> GetGuiInterface() const{
		return fGuiInterface;
	}
	
	bool AmIKilled(){
		if(fGuiInterface){
			return fGuiInterface->AmIKilled();
		}
		else return false;
	}
	
	/**
	 * @brief Set the computational mesh of the analysis.
	 **/
	void SetCompMesh(TPZCompMesh * mesh);
	
	/**
	 * @brief Create an empty TPZAnalysis object
	 **/
	TPZAnalysis();
	
	/**
	 * @brief Destructor: deletes all protected dynamic allocated
	 *objects
	 **/
	virtual ~TPZAnalysis(void);
	/**
	 * @brief Sets the computer connection block number from the graphical
	 *connections block number otimization
	 **/
	void SetBlockNumber();
	
	/**
	 * @brief Returns the dimension of the material which has
	 *the highest dimension
	 **/
	int HighestDimension();
	
	/**
	 * @brief Recompute the node sequence
	 **/
	void Resequence(int firstel = -1);
	
	/**
	 * @brief Assemble the stiffness matrix
	 **/
	virtual  void Assemble();
	
	virtual void AssembleResidual();
	
	/**
	 * @brief Invert the stiffness matrix
	 **/
	virtual void Solve();
	
	/**
	 * @brief Returns the load vector
	 **/
	TPZFMatrix<REAL> &Rhs() { return fRhs;}
	
	
	/**
	 * @brief Returns the solution matrix
	 **/
	TPZFMatrix<REAL> &Solution() { return fSolution;}
	
	/**
	 * @brief Returns the pointer to the computational mesh
	 **/
	TPZCompMesh *Mesh()const { return fCompMesh;}
	/**
	 * @brief Returns a reference to the structural matrix
	 */
	TPZAutoPointer<TPZStructMatrix> StructMatrix() { return fStructMatrix;}
	
	/**
	 * @brief Define the type of preconditioner used
	 *
	 * This method will create the stiffness matrix but without assembling
	 */
	TPZMatrixSolver<REAL> *BuildPreconditioner(EPrecond preconditioner, bool overlap);
	
	void SetTime(REAL time);
	REAL GetTime();
	
private:
	
	/**
	 * @brief Build a sequence solver based on the block graph and its colors
	 */
	TPZMatrixSolver<REAL> *BuildSequenceSolver(TPZVec<int> &graph, TPZVec<int> &graphindex, int neq, int numcolors, TPZVec<int> &colors);
	
	
public:
	void ShowShape( TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames,//1o : TPZConnect* nod,
				   char *plotfile, std::ostream &out=std::cout);
	
	void LoadShape(double dx,double dy, int numelem,TPZConnect* nod);
	
	/** @brief Calls the appropriate sequence of methods to build a
	 * solution or a time stepping sequence
	 */
	virtual void Run(std::ostream &out = std::cout);
	
	virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile);
	
	virtual void CloseGraphMesh();
	
	/**
	 * @brief Defines the postprocessing parameters for the graphical grid
	 **/
	TPZGraphMesh *GraphMesh(int dimension) {
		return fGraphMesh[dimension-1];
	}
	
	virtual void PostProcess(int resolution);
	
	virtual void PostProcess(int resolution, int dimension);
	
	virtual void DefineElementTable(int dimension, TPZVec<int> &GeoElIds, TPZVec<REAL> &points);
	
	virtual void SetTablePostProcessFile(char *filename);
	
	virtual void SetTableVariableNames(int numvar, char **varnames);
	
	virtual void PrePostProcessTable();
	
	virtual void PostProcessTable();
	/**
	 * @brief Load the solution into the computable grid
	 */
	virtual void LoadSolution();
	
	virtual void LoadSolution(const TPZFMatrix<REAL> &sol){
		this->Solution() = sol;
		this->LoadSolution();
	}
	
	void SetExact(void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix<REAL> &deriv));
	
	virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
	
	void PostProcessTable(  TPZFMatrix<REAL> &pos,std::ostream &out= std::cout );
	
	void Print( const std::string &name , std::ostream &out );
	
	TPZMatrixSolver<REAL> & Solver();
	
	void AnimateRun(int num_iter, int steps,
					TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames, const std::string &plotfile);
	
	void SetSolver(TPZMatrixSolver<REAL> &solver);
	
	void SetStructuralMatrix(TPZAutoPointer<TPZStructMatrix> strmatrix);
	
	void SetStructuralMatrix(TPZStructMatrix &strmatrix);
	
};


inline void

TPZAnalysis::SetExact(void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix<REAL> &deriv))
{
	fExact=f;
}

inline TPZMatrixSolver<REAL> &

TPZAnalysis::Solver(){
	return (*fSolver);
}

inline void TPZAnalysis::SetTime(REAL time){
	this->fTime = time;
}

inline REAL TPZAnalysis::GetTime(){
	return this->fTime;
}

#endif
