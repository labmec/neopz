/**
 * @file
 * @brief Contains declaration of TPZMesh class which defines a geometrical mesh and contains a corresponding lists of elements, nodes and conditions.
 */

#ifndef PZGEOMESHH
#define PZGEOMESHH


#include "TPZSavable.h"
#include "pzreal.h"
#include "pzeltype.h"
#include "pzgnode.h"
#include "pzadmchunk.h"
#include "tpzautopointer.h"
#include "pzgeoel.h"

#include <iostream>
#include <string>
#include <map>
#include <list>

/** @brief Identifier indicating the no material is associated */
#define GMESHNOMATERIAL -9999

class TPZMaterial;
class TPZGeoNode;
class TPZCompMesh;
class TPZRefPattern;
class TPZRefPatternDataBase;

template<class T>
class TPZVec;

template <class TGeo> class TPZGeoElRefPattern;

/**
 * @brief This class implements a geometric mesh for the pz environment. \ref geometry "Geometry"
 * @ingroup geometry
 */
/**
 * A TPZGeoMesh defines a geometrical mesh and contains a corresponding list of geometrical elements, geometrical nodes,
 * elementwise defined boundary conditions and nodal boundary conditions. \n
 * Various methods are defined to add, delete or loop over the items which are contained within the TPZGeoMesh. \n
 * Other auxiliary data structures help in the construction of the mesh
 */

class  TPZGeoMesh : public TPZSavable {
	
protected:
	/** @brief TPZGeoMesh name for model identification */
	std::string fName;
	
	/** @brief Computational mesh associated */
	TPZCompMesh 	*fReference;
	
	/** @brief List of pointers to finite elements */
	TPZAdmChunkVector<TPZGeoEl *> fElementVec;
	
	/** @brief List of nodes */
	TPZAdmChunkVector<TPZGeoNode> fNodeVec;
	
	/** @brief Maximum id used by all nodes of this mesh */
	int64_t fNodeMaxId;
	
	/** @brief Maximum id used by all elements of this mesh */
	int64_t fElementMaxId;
    
    /** @brief dimension of the geometric domain */
    int fDim;
	
	typedef std::map<std::pair<int,int>, int> InterfaceMaterialsMap;
	
	/**
	 * @brief Datastructure which indicates the index of the interface material which needs to be created between two different materials 
	 * @see AddInterfaceMaterial 
	 */
	InterfaceMaterialsMap fInterfaceMaterials;
	
public:
	/** @brief Constructors and destructor */
	TPZGeoMesh();
	
	/** @brief Copy constructor */
	TPZGeoMesh(const TPZGeoMesh &cp);
	
	/** @brief Operator of copy */
	TPZGeoMesh & operator= (const TPZGeoMesh &cp );
	
	/** @brief Destructor */
	virtual ~TPZGeoMesh();
	
	/** @brief Deletes all items in the TPZGeoMesh */
	void CleanUp();
	
	/** @brief Reset all connectivities */
	void ResetConnectivities();
	
	public:
virtual int ClassId() const;

	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Indicates that a node with id was created */
	void SetNodeIdUsed(int64_t id) { fNodeMaxId = (id > fNodeMaxId) ? id : fNodeMaxId; }
	
	/** @brief Indicates that an element with id was created */
	void SetElementIdUsed(int64_t id) { fElementMaxId = (id > fElementMaxId) ? id : fElementMaxId; }
	
	/** @brief Returns ++fNodeMaxId */
	int64_t CreateUniqueNodeId() { return ++fNodeMaxId; }
	
	/** @brief Returns ++fElementMaxId */
	int64_t CreateUniqueElementId() { return ++fElementMaxId; }
	
	/** @brief Used in patch meshes */
	void SetMaxNodeId(int64_t id) { fNodeMaxId = (id > fNodeMaxId) ? id : fNodeMaxId; }
	
	/** @brief Used in patch meshes */
	void SetMaxElementId(int64_t id) { fElementMaxId = (id > fElementMaxId) ? id : fElementMaxId; }
	
	/** @brief Number of nodes of the mesh */
	int64_t NNodes() const {return fNodeVec.NElements();}
	
	/** @brief Number of elements of the mesh */
	int64_t NElements() const {return fElementVec.NElements();}
	
	int64_t ReallyNEl() {return (fElementVec.NElements() - fElementVec.NFreeElements()) ; }
	
	void SetName(const std::string &nm);
	
	std::string &Name() { return fName; }
	
	/** @brief Methods for handling pzlists */
	TPZAdmChunkVector<TPZGeoEl *> &ElementVec() { return fElementVec; }
    TPZGeoEl * Element(int64_t iel) { return fElementVec[iel]; }
	TPZAdmChunkVector<TPZGeoNode> &NodeVec() { return fNodeVec; }
	const TPZAdmChunkVector<TPZGeoEl *> &ElementVec() const { return fElementVec; }
	const TPZAdmChunkVector<TPZGeoNode> &NodeVec() const { return fNodeVec; }

    /// Compute the area of the domain
    REAL Area();
    
    /// Compute the area of the domain
    REAL Area(int matid);
    
    /// Compute the area of the domain
    REAL Area(std::set<int> &matids);
    
	/** @brief Resets all load references in elements and nodes */
	void ResetReference();
	
	/** @brief Restore all reference in elements from computational mesh criated from current
     geometrical mesh previously */
	void RestoreReference(TPZCompMesh *cmesh);
	
	/** @brief Sets the reference of the geometric grid to ref */
	void SetReference(TPZCompMesh *ref)
	{
		fReference = ref;
	}
	
	/** @brief Returns the currently loaded computational grid */
	TPZCompMesh *Reference() const {return fReference;}
	
	/** @brief Print the information of the grid to an ostream */
	virtual void Print(std::ostream & out = std::cout);
	virtual void PrintTopologicalInfo(std::ostream & out = std::cout);
	
    /** @brief Returns the nearest node to the coordinate. This method is VERY INEFFICIENT */
	TPZGeoNode* FindNode(TPZVec<REAL> &co);
    
    TPZGeoNode* FindNode(TPZVec<REAL> &co, int &nodeFoundIndex);
    
    /** by Phil */
    /** @brief Returns the element that contains the given point x and it respective point in parametric domain qsi */
    TPZGeoEl * FindElement(TPZVec<REAL> &x, TPZVec<REAL> & qsi, int64_t & InitialElIndex, int targetDim) const;
    
    /** by Caju */
    /** @brief Returns the element that contains the given point x and it respective point in parametric domain qsi */
    TPZGeoEl * FindElementCaju(TPZVec<REAL> &x, TPZVec<REAL> & qsi, int64_t & InitialElIndex, int targetDim);
    
    /** @brief find an element/parameter close to the point */
    TPZGeoEl *FindApproxElement(TPZVec<REAL> &x, TPZVec<REAL> & qsi, int64_t & InitialElIndex, int targetDim) const;
    
    /** by Caju 2013 */
    /** @brief Returns the subelement that contains the given point x and it respective point in parametric domain qsi */
    TPZGeoEl * FindSubElement(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> & qsi, int64_t & InitialElIndex) const;
	
    /** by Philippe 2013 */
    /** @brief Returns the element that is close to the given point x */
    TPZGeoEl * FindCloseElement(TPZVec<REAL> &x, int64_t & InitialElIndex, int targetDim) const;
	
	/** @brief Alternative method for computing the connectivity */
	void BuildConnectivityOld();
	
	/** @brief Build the connectivity of the grid */
	void BuildConnectivity();
	
	/** @brief Fills the nodep vector with pointers to the nodes identified by their indexes */
	void GetNodePtr(TPZVec<int64_t> &nos,TPZVec<TPZGeoNode *> &nodep);
	
	/**
	 * @brief GetBoundaryElements returns all elements beweeen NodFrom and NodTo counterclock wise this
	 * method uses the connectivity of the elements BuildConnectivity should be called to initialize
     */
	/**
	 * The connectivity information this method will only work for grid with 2-D topology the current 
	 * version will only work for a grid with only one level
	 */
	void GetBoundaryElements(int64_t IndexNodeFrom,int64_t IndexNodeTo,TPZStack<TPZGeoEl *> &ElementVec,TPZStack<int> &Sides);
	
    /** @brief Find the element with ids elid */
	TPZGeoEl * FindElement(int64_t elid);
	
	/** @brief Returns the index of the given element into the fElementVec */
	/** @since 2002-05-02 (Cesar) */
	int64_t ElementIndex(TPZGeoEl * gel);
	
	/** @brief Returns the index of the given node into the fNodeVec */
	/** @since 2002-05-02 (Cesar) */
	int64_t NodeIndex(TPZGeoNode * nod);
	
	/**
	 * @brief Generic method for creating a geometric element. Putting this method centrally facilitates
	 * the modification of the element type all through the code
	 * @param type element topology
	 * @param matid material id
	 * @param cornerindexes indexes of the corner nodes of the element
	 * @param index index of the element in the vector of element pointers
	 * @param reftype defines the type of refinement : 0 -> uniform 1-> refinement pattern
	 */
	virtual  TPZGeoEl *CreateGeoElement(MElementType type,TPZVec<int64_t> &cornerindexes,int matid,int64_t &index, int reftype = 1);
	
	/** @brief Creates a geometric element in same fashion of CreateGeoElement but here the elements are blend, as Caju master thesis */
	virtual TPZGeoEl *CreateGeoBlendElement(MElementType type, TPZVec<int64_t>& nodeindexes, int matid, int64_t& index);
	
	/**
	 * @brief Centralized method to delete elements
	 * @param gel pointer to the element to be deleted
	 * @param index index of the element
	 */
	void DeleteElement(TPZGeoEl *gel, int64_t index = -1);
	
	/**
	 * @brief Add an interface material associated to left and right element materials.
	 * @since Feb 05, 2004
	 */
	/**
	 * If std::pair<left, right> already exist, nothing is done and method returns 0.
	 * If material is inserted in geomesh method returns 1.
	 */
	int AddInterfaceMaterial(int leftmaterial, int rightmaterial, int interfacematerial);
	
	/**
	 * @brief Returns the interface material associated to left and right element materials.
     * If no interface material is found GMESHNOMATERIAL is returned
	 * @since Feb 05, 2004
	 */
	int InterfaceMaterial(int leftmaterial, int rightmaterial);
	
	/**
	 * @brief Delete all interface materials in map.
	 * @since Feb 05, 2004
	 */
	void ClearInterfaceMaterialsMap();
    
    /**
	 * @brief Set Dimension.
	 * @since April 17, 2004
	 */
    void SetDimension(int dim){fDim = dim;}
    
    /**
	 * @brief Get Dimension.
	 * @since April 17, 2004
	 */
    int Dimension(){return fDim;}
	
private:
	
	/** @brief Find all elements in elmap or neighbour of elements in elmap which contain a node */
	void BuildElementsAroundNode(int64_t currentnode,std::map<int64_t,TPZGeoEl *> &elmap);
	
	/** @brief Method which works only for two dimensional topologies! */
	/** Find, within elmap, the element which has currentnode as its first boundary side node */
 	void FindElement(std::map<int64_t,TPZGeoEl *> &elmap,int64_t currentnode,TPZGeoEl* &candidate,int &candidateside);
};

#endif
