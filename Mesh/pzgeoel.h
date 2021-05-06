/**
 * @file TPZGeoEl
 * @brief Contains declaration of TPZGeoEl class which defines the behaviour of geometric element.
 */

#ifndef GEOELEMHPP
#define GEOELEMHPP

#include <iostream>

#include "TPZSavable.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzeltype.h"
#include "pztrnsform.h"
#include "doxmesh.h"
#include "pzfmatrix.h"
#include "Hash/TPZHash.h"

#include "pzgeoelside.h"
#include "fadType.h"

class TPZGeoNode;
class TPZCompMesh;
class TPZCompEl;
class TPZGeoMesh;
class TPZCompElSide;
class TPZIntPoints;
class TPZRefPattern;
template<class T>
class TPZVec;
template<class T, int N>
class TPZStack;

/**
 * @brief Defines the behaviour of all geometric elements. \ref geometry "Geometry"
 * @ingroup geometry
 * TPZGeoEl is the common denominator for all geometric elements.
 */

class TPZGeoEl : public virtual TPZSavable {
	
protected:
	
	/** @brief Pointer to the mesh to which the element belongs */
	TPZGeoMesh *fMesh;
    
	/** @brief Traditional element number or element id */
	int64_t		fId;
    
	/** @brief Material index*/
	int		fMatId;
	
    /** @brief Reference to the element currently loaded. Pointer is given as this->Mesh()->Reference()->ElementVec()[fReference] */
	TPZCompEl * fReference;
	
	/** @brief Index of the element from which the element is a subelement*/
	int64_t fFatherIndex;
    
	/** @brief Index of the element in the element vector */
	int64_t fIndex;
    
	/** @brief 3x3 unit matrix to be copied to the axes if the geometric element does not have a particular orientation*/
	static TPZFMatrix<REAL> gGlobalAxes;
    
	/** @brief A counter to indicate how many interface elements are pointing to it */
	int fNumInterfaces;
	
public:
	
    /** Returns the eldest ancestor of this geoel */
    TPZGeoEl * EldestAncestor() const;
    
    /** Returns the directions of this geoel */
//    virtual void Directions(int side,TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions, TPZVec<int> &vectorsides)  = 0;
    
    /** Returns the directions of master element */
    virtual void HDivDirectionsMaster(TPZFMatrix<REAL> &directions) = 0;
    
    /** Returns the directions of this geoel */
    virtual void HDivDirections(TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions)  = 0;

    /** Returns the directions of this geoel */
    virtual void HDivDirections(TPZVec<REAL> &pt, TPZFMatrix<Fad<REAL> > &directions)  = 0;

    /** Returns the eldest ancestor of this geoel */
	virtual void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<REAL> &trans) = 0;
	
	/** @brief Returns number of TPZInterfaceElement pointing to this */
	int NumInterfaces(){
		return this->fNumInterfaces;
	}
	
	/** @brief Increments number of TPZInterfaceElement pointing to this */
	int IncrementNumInterfaces(){
		this->fNumInterfaces++;
		return this->fNumInterfaces;
	}
	 
	/** @brief Decrement number of TPZInterfaceElement pointing to this */
	int DecrementNumInterfaces(){
		this->fNumInterfaces--;
		return this->fNumInterfaces;
	}

	/**
	 * This method gets the ith valid permutation of its topology
	 * @param i number of the permutation to get
	 * @param permutation vector contained the permuted sides
	 */
	virtual void GetPermutation(const int& i, TPZVec<int> &permutation) const = 0;

	/**
	 * @brief Creates an integration rule for the topology of the corresponding side
	 * and able to integrate a polynom of order exactly
	 */
	virtual TPZIntPoints * CreateSideIntegrationRule(int side, int order) =0;
	
	/**
	 * @brief Computes the values of the "num" one dimensional shapefunctions
	 * and derivatives at point x, using lagrangian interpolation
	 * @param x point at which the shapefunction is computed
	 * @param num number of shapefunctions
	 * @param phi values of the shapefunctions
	 * @param dphi values of the derivatives of the shapefunctions
	 */
	void Shape1d(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/**
	 * @brief Computes the values of the "num" one dimensional shapefunctions
	 * at point x, using lagrangian interpolation
	 * @param x point at which the shapefunction is computed
	 * @param num number of shapefunctions
	 * @param phi values of the shapefunctions
	 */
	void ShapePhi1d(REAL x,int num,TPZFMatrix<REAL> &phi);
	
	/** 
	 * @brief Constructor
	 * @param id is the number of the element
	 * @param materialindex is the material index
	 * @param mesh is a pointer to the mesh to which the element belongs
	 */
	TPZGeoEl(int64_t id,int materialindex,TPZGeoMesh &mesh);
	/** 
	 * @brief This constructor generates a unique Id
	 * @param materialindex is the material index
	 * @param mesh is a pointer to the mesh to which the element belongs
	 */
	TPZGeoEl(int materialindex,TPZGeoMesh &mesh);
	
	/** 
	 * @brief This constructor generates a unique Id
	 * @param materialindex is the material index
	 * @param mesh is a pointer to the mesh to which the element belongs
	 * @param index index of the new element in the element vector
	 */
	TPZGeoEl(int materialindex,TPZGeoMesh &mesh,int64_t &index);
	
	/** @brief Copy constructor */
	TPZGeoEl(const TPZGeoEl &el) ;
	
	/** @brief Copy constructor with elements in different meshes */
	TPZGeoEl(TPZGeoMesh & DestMesh, const TPZGeoEl &cp);
	
	/** @brief Copy constructor to a patch mesh */
	TPZGeoEl(TPZGeoMesh & DestMesh, const TPZGeoEl &cp, std::map<int64_t,int64_t> &org2clnMap);
	
	TPZGeoEl() : TPZRegisterClassId(&TPZGeoEl::ClassId), fMesh(0), fId(-1), fMatId(0), fReference(NULL), fFatherIndex(-1), fIndex(-1), fNumInterfaces(0) {
	}
	
	virtual void Initialize()
	{
	}
	
        static int StaticClassId();
        
    int ClassId() const override;
        
	void Read(TPZStream &str, void *context) override;
	
	void Write(TPZStream &str, int withclassid) const override;
	
	virtual TPZGeoEl * Clone(TPZGeoMesh &DestMesh) const = 0;
	
    void GetNodeIndices( TPZVec<int64_t> &nodeindices );
    
    void GetNodeIndices( std::set<int64_t> &nodeindices );
    
	/**
	 * @brief Creates a clone of this element into a new patch mesh
	 * @param DestMesh destination patch mesh
	 * @param gl2lcNdIdx map between node indexes in original and patch mesh
	 * @param gl2lcElIdx map between element indexes in original and patch mesh
	 */
	/** 
	 * This method differs from the above by the patch mesh does not include all elements and nodes. \n
	 * Therefore, a map between node indexes in both meshes are required
	 */
	virtual TPZGeoEl * ClonePatchEl(TPZGeoMesh &DestMesh,
									std::map<int64_t,int64_t> &gl2lcNdIdx,
									std::map<int64_t,int64_t> &gl2lcElIdx) const = 0;
	
	/** @brief Destructor*/
	virtual ~TPZGeoEl();
	
	/** @brief It removes the connectivities of the element*/
	void RemoveConnectivities();
	
	/** @name data access methods
	 * @brief Methods which allow to access the internal data structure of the element
	 */
	//@{
	
	/** @brief Returns the mesh to which the element belongs*/
	TPZGeoMesh *Mesh() const { return fMesh;}
	
	/** @brief Returns the Id of the element*/
	int64_t Id() const { return fId; }
	
	/** @brief Returns the number of nodes of the element*/
	virtual int NNodes() const = 0;
	
	/** @brief Returns the number of corner nodes of the element*/
	virtual int NCornerNodes() const = 0;

    /**
     * Get the number of valid permutations among the element nodes
     * @return
     */
    virtual int NPermutations() const = 0;
	
	/** @brief Returns a pointer to the ith node of the element*/
	TPZGeoNode* NodePtr(int i) const;

	/** @brief Returns the ith node of the element*/
	TPZGeoNode& Node(int i) const;
    
	/**
	 * @brief Returns the index of the ith node the index is the location of the node
	 * in the nodevector of the mesh
	 */
	virtual int64_t NodeIndex(int i) const = 0;

	/** @brief Returns the material index of the element*/
	int MaterialId() const { return fMatId; }

	/** @brief Return a pointer to the element referenced by the geometric element*/
	TPZCompEl *Reference() const;
    
    virtual void GetHigherSubElements(TPZVec<TPZGeoEl*> &unrefinedSons);

	/** @brief Returns the element type acording to pzeltype.h */
	virtual MElementType Type() const = 0;
	/** @brief Returns the element type of a side acording to pzeltype.h */
	virtual MElementType Type(int side) const = 0;
	/** @brief Returns the type of the element as a string */
    virtual std::string TypeName() const
	{
		std::cout << "ElementType should never be called\n";
		return "Notype";
	}
	
	virtual bool IsLinearMapping() const
	{
		return IsLinearMapping(NSides()-1);
	}
    
    virtual bool IsLinearMapping(int side) const = 0;
	
	virtual bool  IsGeoBlendEl() const
	{
		return false;
	}
	
	/** @brief Returns if is a TPZGeoElMapped< T > element */
	/** It is necessary due to the lack of dynamic cast for these elements */
	virtual bool IsGeoElMapped() const{
		return false;
	}
	
	/** @brief Returns the number of connectivities of the element*/
	virtual int NSides() const = 0;
	
	/** @brief Returns the number of nodes for a particular side*/
	virtual int NSideNodes(int side) const = 0;
    
    /** @brief Returns the number of sides of a given dimension */
    virtual int NSides(int dim) const = 0;
	
    /** @brief Returns the first side of a given dimension */
    virtual int FirstSide(int dim) const = 0;
	
	/** @brief Returns the pointer to the nodenum node of side*/
	virtual TPZGeoNode *SideNodePtr(int side,int nodenum) const;
	
	/** @brief Returns the midside node index along a side of the element*/
	virtual void MidSideNodeIndex(int side,int64_t &index) const = 0;
	
	/** @brief Returns the midside node indices along a side of the element */
	/**
	 * THIS METHOD SHOULD SUBSTITUTE MidSideNodeIndex in the future as it is ready for Refinement patterns \n
	 * whereas the former is not
	 *///@TODOFran: this method should be removed, as it is misleading. there is only one midsidenode allowed
	virtual void MidSideNodeIndices(int side,TPZVec<int64_t> &indices) const;
	
	/** @brief Returns the index of the nodenum node of side*/
	virtual int64_t SideNodeIndex(int side,int nodenum) const = 0;
	
	/** @brief Returns the local index of a node on a side*/
	virtual int SideNodeLocIndex(int side, int nodenum) const = 0;
	
	/** @brief Computes the permutation for an HDiv side */
	virtual void HDivPermutation(int side, TPZVec<int> &permutegather);
	
	/** @brief Returns 1 if the side has not been defined by buildconnectivity */
	/**
	 * After construction the side is undefined. The buildconnectivity method loops over all elements \n
	 * and tries to identify neighbours along their uninitialized sides
	 */
	virtual int SideIsUndefined(int side) = 0;
	
	/**
	 * @brief Returns the number of subelements of the element independent of the fact \n
	 * whether the element has already been refined or not 
	 */
	virtual int NSubElements() const = 0;

	/** @brief Returns the number of subelements as returned by GetSubElements2(side) */
	virtual int NSideSubElements(int side) const = 0;
	
	/// Computes the normal vectors needed for forming HDiv vector valued shape functions
//    virtual void VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )=0;
	
	/** @brief Returns a pointer to the father*/
	TPZGeoEl *Father() const;
    
	/** @brief Returns a pointer to the higher level father*/
    TPZGeoEl *LowestFather()
	{
		TPZGeoEl * lowFather(this);
        while(lowFather->Father())
        {
            lowFather = lowFather->Father();
        }
        
        return lowFather;
	}
    
	
	int64_t FatherIndex() { return fFatherIndex; }
	
	/// Set connectivity information elements with blend geometric map
	void BuildBlendConnectivity();

    /**
	 * If the element is a TPZGeoBlend element, this method will ensure that if the side side is connected to the
	 * element with index index, its blend connectivity is erased (for instance, this element may have been deleted).
	 * If it is not a TPZGeoBlend element, the method will just return false.
	 * @param side side in which to seek for connectivities
	 * @param index index of the element that will be disconnected from this
	 * @return true if the element is a TPZGeoBlend element.
	 */
    virtual bool ResetBlendConnectivity(const int64_t &side, const int64_t &index) = 0;
	/** @brief TPZGeoBlend need to find a non-linear neighbour */
	/** 
	 * Although TPZGeoElMapped have non-linear mapping
	 * they are less suitable than self non-linear elements
	 */
	void SetNeighbourForBlending(int side);
	
	//@}

	/**
	 * @brief Method which creates a computational boundary condition element based on the current geometric element, \n
	 * a side and a boundary condition number
	 */
	virtual TPZCompEl *CreateBCCompEl(int side, int bc, TPZCompMesh &cmesh);
	
	/** @brief Creates a geometric element according to the type of the father element */
	virtual TPZGeoEl *CreateGeoElement(MElementType type,
                                       TPZVec<int64_t>& nodeindexes,
                                       int matid,
                                       int64_t& index) = 0;
	
	/** @brief Method which creates a geometric element on the side of an existing element */
	virtual TPZGeoEl *CreateBCGeoEl(int side, int bc) = 0;
	
	/** @brief Returns the side number which is connected to the SideNodes
     returns -1 if no side is found*/
	int WhichSide(TPZVec<int64_t> &SideNodeIds);
	
	/** @brief Returns 1 if gel is a neighbour of the element along side*/
	int NeighbourExists(int side,const TPZGeoElSide &gel);
	
	/** @brief Sets the material index of the element*/
	void SetMaterialId(int id) { fMatId = id;}
	
	/** @brief Initializes the node i of the element*/
	virtual void SetNodeIndex(int i,int64_t nodeindex) = 0;
	
	/** @brief Flags the side as defined, this means no neighbouring element was found*/
	virtual void SetSideDefined(int side) = 0;
	
	/** @brief Returns a pointer to the neighbour and the neighbourside along side of the current element*/
	virtual TPZGeoElSide Neighbour(int side) = 0;
	
	/** @brief Returns the neighbour index for a given side*/
	virtual int64_t NeighbourIndex(int side) const = 0;
	
	/** @brief Fill in the data structure for the neighbouring information*/
	virtual void SetNeighbour(int side,const TPZGeoElSide &neighbour) = 0;
	
	/** @brief Print all relevant data of the element to cout*/
	virtual  void Print(std::ostream & out = std::cout);
    
	/**
	 * @brief Prints the coordinates of all nodes (geometric)
	 * @param out Indicates the device where the coordinates will be printed
	 */
	 virtual void PrintTopologicalInfo(std::ostream & out = std::cout);
	
	/** @brief Make the current element reference to the computational element*/
	void SetReference(TPZCompEl *elp);
	
	/** @brief Sets the subelement of index i */
	virtual void SetSubElement(int i, TPZGeoEl* gel) = 0;
	
	/** @brief Initializes the external connectivities of the subelements */
	virtual void SetSubElementConnectivities();
	
	/** @brief Reset the element referenced by the geometric element to NULL */
	void ResetReference() { this->fReference = NULL; }
	
	/** @brief Equivalent to Print */
	/** @note Overload operator << to print geometric element data */
	friend std::ostream& operator<<(std::ostream &s,TPZGeoEl &el);
	
	/** @brief Divides the element and puts the resulting elements in the vector */
	virtual void Divide(TPZVec<TPZGeoEl *> &pv);
	
    /** @brief Generates a random point in the master domain */
    virtual void RandomPoint(TPZVec<REAL> &pt) = 0;

    /** @brief Return 1 if the element has subelements */
	virtual int HasSubElement() const = 0;

	/**
	 * @brief Compute the transformation between the master element space of one side of an element 
	 * to the master element space of a higher dimension side
	 */
	virtual TPZTransform<REAL> SideToSideTransform(int sidefrom,int sideto)= 0;
    
    /// Project the point from one side to another. The dimension of the points needs to be configured properly
    virtual void ProjectPoint(int sidefrom, TPZVec<REAL> &ptin, int sideto, TPZVec<REAL> &ptout)
    {
        TPZTransform<REAL> tr = SideToSideTransform(sidefrom, sideto);
        tr.Apply(ptin, ptout);
    }
	
	/** @brief Compute the projection of the point within the interior of the element to the side of the element */
	TPZTransform<REAL> Projection(int side);
	
	void SetIndex(int64_t index)
	{
		fIndex = index;
	}
    
    void SetId(int64_t elId)
    {
        fId = elId;
    }
	
	/** @brief Get the transform id the face to face*/
	int GetTransformId2dQ(TPZVec<int> &idfrom,TPZVec<int> &idto);
	
	/** @brief Get the transform id the face to face*/
	int GetTransformId2dT(TPZVec<int> &idfrom,TPZVec<int> &idto);
	
	virtual	TPZTransform<REAL> GetTransform(int side,int son) = 0;

	/** @brief Sets the father element*/
	void SetFather(TPZGeoEl *father)
	{
	    if(!father) fFatherIndex = -1;
		else  fFatherIndex = father->Index();
	}

	/** @brief Sets the father element index
	 * This method is not called SetFather in order to avoid implicit conversions from nullptr to int
	 * */
	virtual void SetFatherIndex(int64_t fatherindex)
	{
		fFatherIndex = fatherindex;
	}
    
    /// return true is gel is an ancestor of the current element
    bool IsSibling(TPZGeoEl *gel)
    {
        if (!gel || fMesh != gel->fMesh) {
            return false;
        }
        TPZGeoEl *father = Father();
        if (father == gel) {
            return true;
        }
        if (father) {
            return father->IsSibling(gel);
        }
        else
        {
            return false;
        }
    }
	
	/** @brief Returns a pointer to the subelement is*/
	virtual TPZGeoEl *SubElement(int is) const = 0;
	
	/** @brief Reset all subelements to NULL*/
	virtual void ResetSubElements()=0;

	/** @brief Returns the number of ancestors */
	int Level();
	
	/** @brief Return the dimension of side*/
	virtual int SideDimension(int side) const = 0;
	
	/** @brief Returns the dimension of the element*/
	virtual int Dimension() const =0;
	
	/** */
	virtual void AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides) = 0;
	virtual void LowerDimensionSides(int side,TPZStack<int> &smallsides) const = 0;
	
    /** @brief Compute a decomposition of the gradient of the mapping function, as a rotation matrix (Jacobian) and orthonormal basis (axes)  */
	void Jacobian(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
    
    /** @brief Compute a decomposition of the gradient of the mapping function, as a rotation matrix (Jacobian) and orthonormal basis (axes)  */
    void JacobianXYZ(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;

    /** @brief Compute a QR facotrization of the gradient of the mapping function, Q = Jacobian and R = axes  */
	static void Jacobian(const TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
    
    /** @brief Compute Jacobian matrix for afine mappings */    
	static void JacobianXYZ(const TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axesXYZ,REAL &detjac,TPZFMatrix<REAL> &jacinv);
    
    /** @brief Return the coordinate in real space of the point coordinate in the master element space*/
    virtual void X(TPZVec<REAL> &qsi,TPZVec<REAL> &result) const = 0;
    
    /** @brief Return the gradient of the transformation at the given coordinate */
    virtual void GradX(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &gradx) const = 0;
    
    /** @brief Return the coordinate in real space of the point coordinate in the master element space*/
    virtual void X(TPZVec<Fad<REAL> > &qsi,TPZVec<Fad<REAL> > &result) const = 0;
    
    /** @brief Return the gradient of the transformation at the given coordinate */
    virtual void GradX(TPZVec<Fad<REAL> > &qsi, TPZFMatrix<Fad<REAL> > &gradx) const = 0;

//	void ComputeNormals(TPZMatrix<REAL> &normal);
	
	/** @brief To test continuity */
	int ElementExists(TPZGeoEl *elem,int64_t id);
	
	/**
	 * @name reftopology
	 * @brief Methods which will implement the declaration of a refinemnt topology
	 * @{
	 */
	
	/** @brief Returns the father/side of the father which contains the side of the sub element */
	virtual TPZGeoElSide Father2(int side) const;
	
	virtual int FatherSide(int side, int son);
	
	/**
	 * @brief Returns the transformation which maps the parameter side of the element/side \n
	 * into the parameter space of the father element/side
	 **/
	virtual TPZTransform<REAL> BuildTransform2(int side, TPZGeoEl *father, TPZTransform<REAL> &t);
	
	
	/**
	 * @brief Returns the side number which is connected to the point pt
	 * @param pt coordinates of the point in parameter space
	 * @return lowest dimension side which contains the point, -1 if no side is found
	 */
	int WhichSide(TPZVec<REAL> &pt);
	
	/** @brief It returns the coordinates from the center of the side of the element in the element coordinate space */
	virtual void CenterPoint(int side, TPZVec<REAL> &masscent) const = 0;
	
	/**
	 * @brief This method will return a partition of the side of the current element \n
	 * as the union of sub elements/side which are put in the stack
	 **/
	virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const;
    
    /**
     * @brief [deprecated] use YoungestChildren
     */
    virtual void GetAllSiblings(TPZStack<TPZGeoEl*> &unrefinedSons);
    /**
     * @brief This method will return all children at the bottom of the refinement tree of the element. i.e. all children that have no subelements
     */
    virtual void YoungestChildren(TPZStack<TPZGeoEl*> &unrefinedSons);
	
	/**
	 * @brief This method will return a partition of the side of the current element \n
	 * as the union of sub elements/side which are put in the stack
	 */
	/** Only element/sides of the given dimension are put on the stack */
	void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel, int dimension) const;
	
	/** @brief Returns the son number of the sub element gel*/
	int WhichSubel() const;
	
	/** @brief Checks the structure of the Father() and GetSubelement2() */
	void CheckSubelDataStructure();

    /**
     * @brief Computes nodes coordinates of the element
     * @param cooridnates xyz coorinates of each node that belongs to element
     */
    void NodesCoordinates(TPZFMatrix<REAL > &cooridnates);
    
	/**
	 * @brief Computes the XInverse and returns true if ksi belongs to master element domain
     * @note ComputeXInverse takes ksi as initial value, so, its recommended that user initialize it
     */
	bool ComputeXInverse(TPZVec<REAL> &XD, TPZVec<REAL> &ksi, REAL Tol);

	/**
	 * @brief Compute the map of a paramenter point in the subelement to a parameter point in the super element
	 * @param ancestor: ancestor element of subelement
	 * @param ksiSon: paramenter point in the subelement
	 * @param ksiAncestor: receives the paramenter point in the ancestor element
	 */
	void TransformSonToFather(TPZGeoEl *ancestor, TPZVec<REAL> &ksiSon, TPZVec<REAL> &ksiAncestor);

	TPZTransform<REAL> ComputeParamTrans(TPZGeoEl *fat,int fatside, int sideson);

	/** @brief Return the volume of the element*/
	REAL Volume();
	
	/** @brief Volume of the master element*/
	virtual REAL RefElVolume() = 0;
	
	/** @brief Returns the area from the face*/
	virtual REAL SideArea(int side);
	
	/** @brief Returns the area from a quadrilateral face*/
	static  REAL QuadArea(TPZVec<TPZGeoNode *> &nodes);
	
	/** @brief Returns the area from the triangular face*/
	static REAL TriangleArea(TPZVec<TPZGeoNode *> &nodes);
	
	virtual REAL ElementRadius();
	
	static REAL Distance(TPZVec<REAL> &centel,TPZVec<REAL> &centface);
	
	/**
	 * @brief Computes the set of normals for defining HDiv approximation spaces
	 * @param normals normal associated with each side
	 * @param vectorsides side associated with each normal vector
	 * @note The normal vectors are initially ordered according to the return of LowerDimensionSides 
	 * and then permuted according to the node id's
	 */
	/** This method will accumulate the normals for all the sides */
//    void ComputeNormals(TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides);
    
    /// Computation of NormalVectors for curvilinear elements
//    void ComputeNormalsDG(TPZVec<REAL> &pt, TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides);
    
    
	virtual REAL CharacteristicSize();
	virtual REAL SmallerEdge();
	/**
	 * @brief Compute the set of normals along a side for defining HDiv approximation spaces
	 * @param side 
	 * @param normals normal associated with each side
	 * @param vectorsides side associated with each normal vector
	 * the normal vectors are initially ordered \n according to the return of LowerDimensionSides
	 * and then permuted according to the node id's
	 */
//    void ComputeNormals(int side, TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides);
//    void ComputeNormalsDG(int side, TPZVec<REAL> &pt, TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides);
	/**
	 * @brief Compute the permutation needed to order the normal vectors in a consistent way
	 * \f$ normal(indexfrom[i]) = normal(i) \f$
	 */
	/** This permutation needs to be applied to the shape functions */
	void ComputePermutationNormals(int side, TPZVec<int> &indexfrom);
	
	/** @brief Determine the orientation of the normal vector comparing the ids of the neighbouring elements */
	int NormalOrientation(int side);
	
	/** @brief Defines the refinement pattern. It's used only in TPZGeoElRefPattern objects. */
	virtual void SetRefPattern(TPZAutoPointer<TPZRefPattern> );
	
	/** @brief Returns the refinement pattern associated with the element */
	virtual TPZAutoPointer<TPZRefPattern> GetRefPattern() const;
	
	/**
	 * @brief Verify coordinate of element nodes checking if they are coincident to the X mapping 
	 * of the corner nodes of parametric elements
	 * @return true if everything OK else false
	 */
	bool VerifyNodeCoordinates(REAL tol = 1e-1);

	/** @brief Verifies if the parametric point pt is in the element parametric domain */
	virtual bool IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol = 1.e-2) = 0;
	
	/**
	 * @brief Ortogonal projection from given qsi to a qsiInDomain (all in the element parametric domain)
	 * @return Returns the side where the point was projected.
     * @note If ortogonal projection to any sides of element results in a qsi outside domain,
     *       this method will returns the nearest node side.
	 * @note Observe that if the point is already in the parametric domain, the method will return \f$ NSides() - 1 \f$
	 */
	virtual int ProjectInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain) = 0;
    
    /**
	 * @brief Projection from given qsi to a qsiInDomain (in the element boundary) using bissection method from given qsi to element center.
	 * @return Returns the side where the point was projected.
	 * @note Observe that if the point is already in the parametric domain, the method will return \f$ NSides() - 1 \f$
	 */
	virtual int ProjectBissectionInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain) = 0;
	
	 /** @brief Returns the index of the element within the element vector of the mesh */
    int64_t Index() const
    {
        return fIndex;
    }
	
private:
	/** @brief To be used after the buid connectivity. If some neighbour isn't initialized */
	/** It will be initialized as the TPZGeoElSide (this, thisside) */
	void InitializeNeighbours();
};


inline void TPZGeoEl::Divide(TPZVec<TPZGeoEl *> &) {
	PZError << "TPZGeoEl::Divide is called.\n";
}

inline void TPZGeoEl::SetReference(TPZCompEl * elp){
	this->fReference = elp;
}

inline TPZCompEl *TPZGeoEl::Reference() const {
	return this->fReference;
}

inline TPZGeoElSide TPZGeoElSide::Neighbour() const {
    return fGeoEl ? fGeoEl->Neighbour(fSide) : TPZGeoElSide();
}



#endif
