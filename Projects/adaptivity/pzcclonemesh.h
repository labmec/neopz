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
class TPZTransfer;
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
    REAL fhError;	//??
    REAL fError;  //??
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
   * return the element index in reference mesh
   */
  int GetOriginalElementIndex(int elindex);


  TPZInterpolatedElement *GetOriginalElement(TPZCompEl *el);

  /**
   * Returns the Reference Mesh
   */
  TPZCompMesh *GetReferenceMesh(){return fCloneReference;}

  /**
   * Evaluates the mesh error as being the difference between reference solution and
   * a solution gived by a uniformly refined mesh
   */
  void MeshError(TPZCompMesh *fine, TPZVec<REAL> &ervec, void(*f)(TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),TPZVec<REAL> &truervec);

  /**
   * Returns the uniformly hp refined mesh base on this mesh
   */
  TPZCompMesh * UniformlyRefineMesh();


  REAL ElementError(TPZInterpolatedElement *fine, 
		    TPZInterpolatedElement *coarse, 
       		    TPZTransform &tr,
		    void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),
		    REAL &truerror);

  /**
   * Returns hp pattern of reference elements
   * @param minerror Minimum error for the element be analysed
   * @param error Vector containing all elements error
   * @param finee Uniformly refined mesh
   * @param gelstack Geometric elements stack. This stack will include the h / hp refined elements
   * @param porder p order of the elements contained in gelstack
   */
  void ApplyRefPattern(REAL minerror, TPZVec<REAL> &error, TPZCompMesh *finee, 
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
  int IsFather(TPZGeoEl *el);


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
  void LoadSolution(TPZFMatrix &sol);


  void Print(std::ostream &out);
};

#endif

