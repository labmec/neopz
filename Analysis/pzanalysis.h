#ifndef ANALYSISH
#define ANALYSISH

class TPZGeoMesh;
class TPZCompMesh;
class TPZMatrix;
class TPZBlock;
class TPZConnect;
class TPZSolver;
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

/**
 * @brief Class TPZAnalysis implements analysis procedures
 * @ingroup analysis
 */
class TPZStructMatrix;

template<class T, int N> class TPZStack;

/// class which implements the sequence of actions to perform a finite element analysis
/**
This class will renumerate the nodes upon construction
*/
class TPZAnalysis {

public:

/// Precondtioners which can be created by objects of this class
 enum EPrecond { EJacobi, EBlockJacobi, EElement, ENodeCentered };


 protected:
  /**
   * Geometric Mesh
   */
  TPZGeoMesh *fGeoMesh;
  /**
   * Computational mesh
   */
  TPZCompMesh *fCompMesh;
  /**
   * Graphical mesh
   */
  TPZGraphMesh *fGraphMesh[3];
  /**
   * Load vector ???
   */
  TPZFMatrix fRhs;
  /**
   * Solution vector
   */
  TPZFMatrix fSolution;
  /**
   * Type of solver to be applied
   */
  TPZMatrixSolver *fSolver;
  /**
   * Scalar variables names - to post process
   */
  TPZVec<std::string> fScalarNames[3];
  /**
   * Vector variables names - to post process
   */
  TPZVec<std::string> fVectorNames[3];
  /**
   * Step ???
   */
  int fStep;
  /**
   * Time variable which is used in dx output
   */
  REAL fTime;

  /**
   * Structural matrix
   */
  TPZAutoPointer<TPZStructMatrix>  fStructMatrix;
  
  /**
   * Renumbering scheme
   */
	TPZAutoPointer<TPZRenumbering> fRenumber;

	/** Pointer for gui interface object */
	TPZAutoPointer<TPZGuiInterface> fGuiInterface;

  /// datastructure which defines postprocessing for one dimensional meshes
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
     * Exact solution ??? not implemented???
     */
    void (*fExact)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix &deriv);

  /**
   *Create an TPZAnalysis object from one mesh pointer
   */
  TPZAnalysis(TPZCompMesh *mesh,std::ostream &out = std::cout);

  /**
   *Create an TPZAnalysis object from one mesh pointer
   */
	TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh,std::ostream &out = std::cout);

	/** Defines gui interface object */
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
	 * Set the computational mesh of the analysis.
	 **/
  void SetCompMesh(TPZCompMesh * mesh);

	/**
   *Create an empty TPZAnalysis object
   **/
  TPZAnalysis();

  /**
   *Destructor: deletes all protected dynamic allocated
   *objects
   **/
  virtual ~TPZAnalysis(void);
  /**
   *Set the computer connection block number from the graphical
   *connections block number otimization
   **/
  void SetBlockNumber();

  /**
   *Returns the dimension of the material which has
   *the highest dimension
   **/
  int HighestDimension();

  /**
   *Recompute the node sequence
   **/
  void Resequence(int firstel = -1);

  /**
   *Assemble the stiffness matrix
   **/
  virtual  void Assemble();
  
  virtual void AssembleResidual();

  /**
   *Invert the stiffness matrix
   **/
  virtual void Solve();

  /**
   *Returns the load vector
   **/
  TPZFMatrix &Rhs() { return fRhs;}


  /**
   *Returns the solution matrix
   **/
  TPZFMatrix &Solution() { return fSolution;}

  /**
   *Returns the pointer to the computational mesh
   **/
  TPZCompMesh *Mesh()const { return fCompMesh;}
  /**
  * Returns a reference to the structural matrix
  */
  TPZAutoPointer<TPZStructMatrix> StructMatrix() { return fStructMatrix;}

  /**
  * Define the type of preconditioner used
  * This method will create the stiffness matrix but without assembling
  */
  TPZMatrixSolver *BuildPreconditioner(EPrecond preconditioner, bool overlap);

  void SetTime(REAL time);
  REAL GetTime();

 private:

 /**
 * Build a sequence solver based on the block graph and its colors
 */
TPZMatrixSolver *BuildSequenceSolver(TPZVec<int> &graph, TPZVec<int> &graphindex, int neq, int numcolors, TPZVec<int> &colors);


public:
  void ShowShape( TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames,//1o : TPZConnect* nod,
		  char *plotfile, std::ostream &out=std::cout);

  void LoadShape(double dx,double dy, int numelem,TPZConnect* nod);

  virtual void Run(std::ostream &out = std::cout);
  // calls the appropriate sequence of methods to build a
  // solution or a time stepping sequence
  virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile);

  virtual void CloseGraphMesh();

  /**
   * defines the postprocessing parameters for the graphical grid
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
   *load the solution into the computable grid
   */
	virtual void LoadSolution();

	virtual void LoadSolution(const TPZFMatrix &sol){
		this->Solution() = sol;
		this->LoadSolution();
	}

  void SetExact(void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix &deriv));

  virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);

  void PostProcessTable(  TPZFMatrix &pos,std::ostream &out= std::cout );

  void Print( const std::string &name , std::ostream &out );

  //misael
  TPZMatrixSolver & Solver();

  void AnimateRun(int num_iter, int steps,
		  TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames, const std::string &plotfile);

  void SetSolver(TPZMatrixSolver &solver);

  void SetStructuralMatrix(TPZAutoPointer<TPZStructMatrix> strmatrix);
  
  void SetStructuralMatrix(TPZStructMatrix &strmatrix);
  
  
  
};


// Inline


inline void

TPZAnalysis::SetExact(void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix &deriv))
{
  fExact=f;
}

inline TPZMatrixSolver &

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

