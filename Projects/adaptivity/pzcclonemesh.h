//HEADER FILE FOR CLASS CLONEMESH

#ifndef PZCCLONEMESHHPP
#define PZCCLONEMESHHPP

#include "pzcmesh.h"
#include "pzonedref.h"
#include <map>

class TPZCompEl;
class TPZGeoEl;
struct TPZCompElBC;
class TPZConnect;
struct TPZConnectBC;
class TPZBndCond;
class TPZMaterial;
class TPZGeoMesh;
class TPZCoSys;
class TPZInterpolatedElement;
class TPZTransform;
class TPZGeoCloneMesh;

	/**
	 * Class TPZCompCloneMesh implemments computational mesh
	 * @ingroup CompMesh
	 */
//-----------------------------------------------------------------------------
// class structure
//-----------------------------------------------------------------------------
class TPZCompCloneMesh : public TPZCompMesh {
  
 protected:
  /**
   * Computational mesh to which this grid is cloned
   */
  TPZCompMesh * fCloneReference;
  
  /**
   * Maps connect index from original mesh to cloned mesh
   */
    std::map<int,int> fMapConnects;
  
  /**
   * Maps connect index from cloned mesh to original mesh
   */
  TPZStack <int> fOriginalConnects;
  
  //  void CleanUp();
  
 public:
  struct TPZRefPattern {
    int fId[3]; 	//Subelements connectivities ids
    int fp[2];		//subelements p-order refinement
    int fh[2];		//subelements h-order refinement
    REAL fhError;	// Error if h refinement is applied
    REAL fError;  // Smallest error
	  TPZRefPattern(int id1,int id2,int id3,int p1,int p2,int h1,int h2,REAL herror,REAL error) {
		  fId[0] = id1; fId[1] = id2; fId[2] = id3;
		  fp[0] = p1; fp[1] = p2;
		  fh[0] = h1; fh[1] = h2;
		  fhError = herror; fError = error;
	  }
	  TPZRefPattern() {
		  fId[0] = fId[1] = fId[2] = -1;
		  fp[0] = fp[1] = fh[0] = fh[1] = 0;
		  fhError = fError = 0.0;
	  }
  };

  /**
   * Copy constructor
   */
  TPZCompCloneMesh (TPZCompCloneMesh &copy);


  /**
   * Constructor from geometrical mesh
   * @param gr pointer to geometrical clone reference mesh
   * @param cmesh pointer to computational refernce mesh
   */
  TPZCompCloneMesh (TPZGeoCloneMesh * gr, TPZCompMesh * cmesh);

  /**
   * Simple Destructor
   */
  virtual ~TPZCompCloneMesh();


  /**
   * @brief Given the index in the CompClone mesh return the comp element index in the reference mesh
   */
  int GetOriginalElementIndex(int elindex);


    /**
     * @brief Given the pointer to the element in the CompClone mesh return the pointer to the comp element in the reference mesh
     */
  TPZInterpolatedElement *GetOriginalElement(TPZCompEl *el);

  /**
   * @brief Returns the Reference Mesh
   */
  TPZCompMesh *GetReferenceMesh(){return fCloneReference;}

  /**
   * Evaluates the mesh error as being the difference between reference solution and
   * a solution computed on a uniformly refined mesh
   * The mesh error is computed for the current clone mesh
   * the error and true error are accumulated in ervec and truervec
   * ervec and truervec are indexed according to the original computational mesh
   */
  void MeshError(TPZCompMesh *fine, TPZVec<REAL> &ervec, void(*f)(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<REAL> &deriv),TPZVec<REAL> &truervec);

  /**
   * Returns the uniformly hp refined mesh base on this mesh
   */
  TPZCompMesh * UniformlyRefineMesh();


  REAL ElementError(TPZInterpolatedElement *fine, 
		    TPZInterpolatedElement *coarse, 
       		    TPZTransform &tr,
		    void (*f)(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<REAL> &deriv),
		    REAL &truerror);

  /**
   * Returns hp pattern of reference elements
   * @param minerror Minimum error for the element be analysed
   * @param error Vector containing all elements error
   * @param fineMesh Uniformly refined mesh
   * @param gelstack Geometric elements stack. This stack will include the h / hp refined elements
   * @param porder p order of the elements contained in gelstack
   */
  void ApplyRefPattern(REAL minerror, TPZVec<REAL> &error, TPZCompMesh *fineMesh, 
		       TPZStack<TPZGeoEl *> &gelstack, TPZStack<int> &porder);


protected:
  /**
   * Verifies if the connect cnid is already mapped
   */
  int HasConnect(int cnid);

  /**
   * Creates the Dirichlet Boundary Condition along the patch sides
   */
  void CreateCloneBC();

  /**
   * Verifies if the given element is son of the Geometric Reference Element
   * @param el - element to analyse
   */
  int IsSonOfRootElement(TPZGeoEl *el);


  /**
   * Analyse an element and return its best refinement hp / p.
   * @param f One dimensional refined structure.
   * @param cint Element to be analysed
   * @param subels Will return the subelements in case of h refinement
   * @param porders Refinement order for each linear side of the element
   */
  void AnalyseElement( TPZOneDRef &f, TPZInterpolatedElement *cint,
		       TPZStack<TPZGeoEl *> &subels, TPZStack<int> &porders);

  void AdaptElements (TPZVec<TPZGeoEl *> &gelstack,TPZVec<int> &porders);


  void DeduceRefPattern(TPZVec<TPZRefPattern> &refpat,	
					  TPZVec<int> &cornerids,
					  TPZVec<int> &porders,
					  int originalp);

  void Sort(TPZVec<REAL> &vec, TPZVec<int> &perm);

  
public:

  //@}

  /**
   * @name SCIENTIFIC_ROUTINES
   * Scientific manipulating routines
   */

  //@{

  /**
   * Creates the computational elements, and the degree of freedom nodes
   */
  virtual void AutoBuild();

  /**
   * Given the solution of the global system of equations, computes and stores the
   * solution for the restricted nodes
   * @param sol given solution matrix
   */
  void LoadSolution(TPZFMatrix<REAL> &sol);


  void Print(std::ostream &out) const;
};

#endif

