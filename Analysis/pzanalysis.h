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
#include <iostream>

using namespace std;

#include "pzfmatrix.h"
/**
 * Class TPZAnalysis implements analysis procedures
 * @ingroup analysis
 */
class TPZStructMatrix;

template<class T, int N> class TPZStack;

class TPZAnalysis {

public:

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
  TPZVec<char *> fScalarNames[3];
  /**
   * Vector variables names - to post process
   */
  TPZVec<char *> fVectorNames[3];
  /**
   * Step ???
   */
  int fStep;
  /**
   * Time ?? (time step ??)
   */
  REAL fTime;

  /**
   * Structural matrix
   */
  TPZStructMatrix * fStructMatrix;

  struct TTablePostProcess {
    TPZVec<int> fGeoElId;
    TPZVec<TPZCompEl *> fCompElPtr;
    int fDimension;
    TPZVec<REAL> fLocations;
    TPZVec<char *> fVariableNames;
    ostream *fOutfile;
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
  TPZAnalysis(TPZCompMesh *mesh,std::ostream &out = cout);

	/**
	 * Set the computational mesh of the analysis.
	 **/
  void TPZAnalysis::SetCompMesh(TPZCompMesh * mesh);

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
  TPZCompMesh *Mesh() { return fCompMesh;}
  /**
  * Returns a reference to the structural matrix
  */
  TPZStructMatrix &StructMatrix() { return *fStructMatrix;}
  
  /**
  * Define the type of preconditioner used
  * This method will create the stiffness matrix but without assembling
  */
  TPZMatrixSolver *BuildPreconditioner(EPrecond preconditioner, bool overlap);
  
 private:
 
 /**
 * Build a sequence solver based on the block graph and its colors
 */
TPZMatrixSolver *BuildSequenceSolver(TPZVec<int> &graph, TPZVec<int> &graphindex, int neq, int numcolors, TPZVec<int> &colors);
  

public:
  void ShowShape( TPZVec<char *> &scalnames, TPZVec<char *> &vecnames,//1o : TPZConnect* nod,
		  char *plotfile, ostream &out=cout);

  void LoadShape(double dx,double dy, int numelem,TPZConnect* nod);

  virtual void Run(ostream &out = cout);
  // calls the appropriate sequence of methods to build a
  // solution or a time stepping sequence
  virtual void DefineGraphMesh(int dimension, TPZVec<char *> &scalnames, TPZVec<char *> &vecnames, const char *plotfile);
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
  void LoadSolution();

  void SetExact(void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix &deriv));

  void PostProcess(TPZVec<REAL> &loc, ostream &out = cout);

  void PostProcessTable(  TPZFMatrix &pos,ostream &out= cout );

  void Print( char *name , ostream &out );

  //misael
  TPZMatrixSolver & Solver();

  void AnimateRun(int num_iter, int steps,
		  TPZVec<char *> &scalnames, TPZVec<char *> &vecnames, char *plotfile);

  void SetSolver(TPZMatrixSolver &solver);

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

         
#endif

