/**
 * @file
 * @brief HEADER FILE FOR CLASS COMPUTATIONAL CLONE MESH
 */

#ifndef PZCCLONEMESHHPP
#define PZCCLONEMESHHPP

#include "pzcmesh.h"
#include "pzonedref.h"
#include <map>

class TPZCompEl;
class TPZGeoEl;
class TPZConnect;
class TPZBndCond;
class TPZMaterial;
class TPZGeoMesh;
class TPZInterpolatedElement;
template<class T>
class TPZTransform;
class TPZGeoCloneMesh;

/**
 * @brief Class TPZCompCloneMesh implements cloned computational mesh
 * @ingroup CompMesh
 */
class TPZCompCloneMesh : public TPZCompMesh {
    
protected:
    /** Computational mesh to which this grid is cloned */
    TPZCompMesh * fCloneReference;
    
    /** Maps connect index from original mesh to cloned mesh */
    std::map<int64_t,int64_t> fMapConnects;
    
    /** Maps connect index from cloned mesh to original mesh */
    TPZStack <int64_t> fOriginalConnects;
        
public:
    struct TPZRefPattern {
        int64_t fId[3]; 	//Subelements connectivities ids
        int fp[2];		//subelements p-order refinement
        int fh[2];		//subelements h-order refinement
        REAL fhError;	// Error if h refinement is applied
        REAL fError;    // Smallest error
        TPZRefPattern(int64_t id1,int64_t id2,int64_t id3,int p1,int p2,int h1,int h2,REAL herror,REAL error) {
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
    
    /** Copy constructor */
    TPZCompCloneMesh(const TPZCompCloneMesh &copy);
    TPZCompCloneMesh();
    
    /**
     * Constructor from geometrical mesh
     * @param gr pointer to geometrical clone reference mesh
     * @param cmesh pointer to computational refernce mesh
     */
    TPZCompCloneMesh(TPZGeoCloneMesh * gr,TPZCompMesh * cmesh);
    
    /** Simple Destructor */
    virtual ~TPZCompCloneMesh();
    
    /**
     * @brief Given the index in the CompClone mesh return the comp element index in the reference mesh
     */
    int64_t GetOriginalElementIndex(int64_t elindex);
    
    /**
     * @brief Given the pointer to the element in the CompClone mesh return the pointer to the comp element in the reference mesh
     */
    TPZInterpolatedElement *GetOriginalElement(TPZCompEl *el);
    
    /** @brief Returns the Reference Mesh */
    TPZCompMesh *GetReferenceMesh() { return fCloneReference; }
    
    /**
     * @brief Evaluates the mesh error as being the difference between reference solution and
     * a solution computed on a uniformly refined mesh
     * The mesh error is computed for the current clone mesh
     * the error and true error are accumulated in ervec and truervec
     * ervec and truervec are indexed according to the original computational mesh
     */
    void MeshError(TPZCompMesh *fine, TPZVec<REAL> &ervec, void(*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),
                   TPZVec<REAL> &truervec, std::ofstream &out);
    
    /**
     * @brief Returns the uniformly hp refined mesh base on this mesh
     * @param maxp Maxime order to using when a element to be divided.
     */
    TPZCompMesh * UniformlyRefineMesh(int maxp,int print=0);
    
    
    REAL ElementError(TPZInterpolatedElement *fine,
                      TPZInterpolatedElement *coarse,
                      TPZTransform<> &tr,
                      void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),
                      REAL &truerror);
    
    /**
     * @brief Returns hp pattern of reference elements
     * @param minerror Minimum error for the element be analysed
     * @param error Vector containing all elements error
     * @param fineMesh Uniformly refined mesh
     * @param gelstack Geometric elements stack. This stack will include the h / hp refined elements
     * @param porder p order of the elements contained in gelstack
     */
    void ApplyRefPattern(REAL minerror, TPZVec<REAL> &error, TPZCompMesh *fineMesh,
                         TPZStack<TPZGeoEl *> &gelstack, TPZStack<int> &porder);
    
	/** @brief Returns the unique identifier for reading/writing objects to streams */
public:
virtual int ClassId() const;

	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
    
    /** @brief To clone this mesh */
   // TPZCompCloneMesh* Clone() const;
    
protected:
    
    /** @brief Verifies if the connect cnid is already mapped */
    int HasConnect(int64_t cnid);
    
    /** @brief Creates the Dirichlet Boundary Condition along the patch sides */
    void CreateCloneBC();
    
	/** @brief Clone this mesh */
	TPZCompCloneMesh * Clone() const;

    /**
     * @brief Verifies if the given element is son of the Geometric Reference Element
     * @param el - element to analyse
     */
    int IsSonOfRootElement(TPZGeoEl *el);
    
    
    /**
     * @brief Analyse an element and return its best refinement hp / p.
     * @param f One dimensional refined structure.
     * @param cint Element to be analysed
     * @param subels Will return the subelements in case of h refinement
     * @param porders Refinement order for each linear side of the element
     */
    void AnalyseElement( TPZOneDRef &f, TPZInterpolatedElement *cint,
                        TPZStack<TPZGeoEl *> &subels, TPZStack<int> &porders);
    
    void AdaptElements (TPZVec<TPZGeoEl *> &gelstack,TPZVec<int> &porders);
    
    
    void DeduceRefPattern(TPZVec<TPZRefPattern> &refpat,
                          TPZVec<int64_t> &cornerids,
                          TPZVec<int> &porders,
                          int originalp);
    
    void Sort(TPZVec<REAL> &vec, TPZVec<int64_t> &perm);
    
    
public:
    
    //@{
    
    /**
     * @name SCIENTIFIC_ROUTINES
     * Scientific manipulating routines
     */
    
    //@}
    
    /** Creates the computational elements, and the degree of freedom nodes and copy solution from original computational mesh */
    virtual void AutoBuild();
    
    /**
     * @brief Given the solution of the global system of equations, computes and stores the
     * solution for the restricted nodes
     * @param sol given solution matrix
     */
    void LoadSolution(TPZFMatrix<REAL> &sol);
    
    
    void Print(std::ostream &out) const;

};

//#ifndef BORLAND
//template class TPZRestoreClass<TPZCompCloneMesh>;
//#endif

#endif
