/**
 * @file
 * @brief Contains declaration of TPZGeoElSide class which represents an element and its side, and TPZGeoElSideIndex class which represents an TPZGeoElSide index.
 */
//$Id: pzgeoelside.h,v 1.28 2011-05-15 21:09:19 phil Exp $

#ifndef PZGEOELSIDEH
#define PZGEOELSIDEH

/*******       TPZGeoElSide       *******/

class TPZGeoEl;
class TPZTransform;
class TPZCompElSide;
template<class TVar>
class TPZFMatrix;

#include "pzvec.h"
#include "pzstack.h"
#include "pzgmesh.h"
#include <set>

class TPZGeoElSide;

/**
 * @ingroup geometry
 * @brief Utility class which represents an element index with its side. \ref geometry "Geometry"
 */
class TPZGeoElSideIndex{
private:
	int fGeoElIndex;
	int fSide;
	
public:
	/** @brief Destructor. */
	~TPZGeoElSideIndex();
	/** @brief Simple constructor. */
	TPZGeoElSideIndex();
	/** @brief Constructor with geometric element referenced and corresponding side. */
	TPZGeoElSideIndex(TPZGeoEl *gel,int side);
	
	TPZGeoElSideIndex(int gelindex,int side);
	
	TPZGeoElSideIndex(const TPZGeoElSide &side);
	
	TPZGeoElSideIndex(const TPZGeoElSideIndex &cp);
	/** @brief To clone current object */
	TPZGeoElSideIndex * Clone();
	/** @brief Redefines operator = attribuition to TPZGeoElSideIndex object. */
	TPZGeoElSideIndex &operator= (const TPZGeoElSideIndex &A );
	
	int Side() const;
	
	void SetSide(int side);
	
	TPZGeoEl *Element(const TPZGeoMesh *mesh) const;
	
	void SetElement(TPZGeoEl* geoel);
	
	int ElementIndex() const;
	
	void SetElementIndex(int i);
	
	void Read(TPZStream &buf);
	void Write(TPZStream &buf);
};

/**
 * @brief Utility class which represents an element with its side. \ref geometry Geometry
 * @ingroup geometry
 */
/**
 * This class is often used to manipulate neighbouring information between elements
 */
class TPZGeoElSide {
	
	TPZGeoEl *fGeoEl;
	int fSide;
public:
	int NSubElements();
	void GetSubElements2(TPZStack<TPZGeoElSide> &subelements);
	TPZGeoElSide StrictFather();
	TPZGeoElSide Father2();
	TPZCompElSide LowerLevelCompElementList2(int onlyinterpolated);
	
	/** @brief Checks whether other is a relative (son or ancestor) of this
	 */
	bool IsRelative(TPZGeoElSide other);
	
	/** @brief Checks whether other is an ancestor of this
	 */
	bool IsAncestor(TPZGeoElSide other);
	
	/** By Caju */
	void X(TPZVec< REAL > &loc, TPZVec< REAL > &result);
	
	/// Jacobian associated with the side of the element
	void Jacobian(TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
    
    /// Area associated with the side
    REAL Area();
	
	int NNeighbours();
	
	/**
	 * @brief Returns the number of neighbours, excluding the given element (thisElem)
	 */
	int NNeighboursButThisElem(TPZGeoEl *thisElem);
	
	/**
	 * Will return all elements of equal or higher level than than the current element
	 * 
	 * All elements/sides have the same dimension
	 */
	//void EqualorHigherCompElementList(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates);
	/**
	 * Will return all elements of equal or higher level than than the current element
	 */
	void EqualorHigherCompElementList2(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates);
	
	TPZGeoElSide(){ fGeoEl = 0; fSide  = -1;}
	//TPZGeoElSide(const TPZGeoElSide &gelside);
	
	TPZGeoElSide(TPZGeoEl *gel,int side){  fGeoEl = gel; fSide = side;}
	
	/**
	 * @brief This constructor set an TPZGeoElSide based in the cornerNodes of an side of gel
	 * 
	 * If the cornerNodes are not consistent, the TPZGeoElSide created is NULL
	 */
	TPZGeoElSide(TPZGeoEl *gel, std::set<int> &sideCornerNodes);
	
	TPZGeoElSide(const TPZGeoElSideIndex &index, const TPZGeoMesh * mesh){
		this->fSide = index.Side();
		this->fGeoEl = index.Element(mesh);
	}//end of method
    
    TPZGeoElSide(int zero) : fGeoEl(0), fSide(-1)
    {
        
    }
	
	TPZGeoEl *Element()const{return fGeoEl;}
	
	void SetElement(TPZGeoEl* geoel){ fGeoEl = geoel; }
	
	int Side() const {return fSide;}
	
	void SetSide(int side) { fSide = side; }
	
	bool IsLinearMapping() const;
	
	int Exists() const {return (fGeoEl != 0 && fSide > -1);}
	
	void CenterPoint(TPZVec<REAL> &center);
    
    /// Returns the number of sides in which the current side can be decomposed
    int NSides();
	
	TPZGeoElSide Neighbour() const;//return neighbour of the side fSide
	
	/**
	 * @brief Returns the set of neighbours which can directly be accessed by the datastructure
	 */
	void AllNeighbours(TPZStack<TPZGeoElSide> &allneigh);
	
	/**
	 * @brief Returns the set of neighbours as computed by the intersection of neighbours along vertices
	 */
	void ComputeNeighbours(TPZStack<TPZGeoElSide> &compneigh);
	
	int Id();
    
    /**
     * @brief the dimension associated with the element/side
     */
	int Dimension() const;
	
	int operator==(const TPZGeoElSide &other) const {
		return fGeoEl == other.fGeoEl && fSide == other.fSide;
	}
	int operator!=(const TPZGeoElSide &other) const {
		return fGeoEl != other.fGeoEl || fSide != other.fSide;
	}
	
	int operator<(const TPZGeoElSide &other) const {
		return (fGeoEl < other.fGeoEl || (fGeoEl == other.fGeoEl && fSide < other.fSide));
	}
	
	int operator>(const TPZGeoElSide &other) const {
		return (fGeoEl > other.fGeoEl || (fGeoEl == other.fGeoEl && fSide > other.fSide));
	}
    
    /**
     * @brief The conversion to bool indicates whether the object has an associated element
     */
    operator bool() const
    {
        return fGeoEl != 0;
    }
	
	//void SideTransform(TPZGeoElSide neighbour,TPZTransform &t);
	/**set the element neighbour along side i equal to el
     the side of the neighbour which is connected to the current element
     is neighbourside*/
	//void SideTransform2(TPZGeoElSide neighbour,TPZTransform &t);
	
	/**
	 @brief Accumulates the transformations from the current element/side to the
	 neighbour/side
	 
	 Third improved version
	 */
	void SideTransform3(TPZGeoElSide neighbour,TPZTransform &t);
	
	void SetConnectivity(const TPZGeoElSide &neighbour) const;
	/**
	 * @brief This method inserts the element/side and all lowerdimension sides into the connectivity loop
	 */
	void InsertConnectivity(TPZGeoElSide &neighbour);
	
	void RemoveConnectivity();
	
	static void BuildConnectivities(TPZVec<TPZGeoElSide> &elvec, TPZVec<TPZGeoElSide> &neighvec);
	
	/** @brief Fill in the data structure for the neighbouring information*/
	void SetNeighbour(const TPZGeoElSide &neighbour) const;
	
	/** Returns a pointer to computational element referenced by a geometric
	 element which is a son along side and has higher level than level
	 
	 if onlyinterpolated = 1 only elements TPZInterpolatedElement will be put on the stack*/
	//  void SmallConnect(int level,TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated);
	
	/**search elements of high dimension connected to the actual element*/
	//TPZCompElSide HigherDimensionConnected(int targetdimension,int onlyinterpolated);
	
	//TPZGeoElSide HigherDimensionSide(int targetdimension);
	
	/**return pointers to the subelements and their sides along side
	 in the order determined by the referencenodes*/
	//void GetSubElement(TPZVec<int> &referencenodes,TPZVec<TPZGeoElSide> &subelements);
	
	/**return all subelements along the current side*/
	//void GetSubElements(TPZVec<TPZGeoElSide> &subelements);
	
	TPZTransform NeighbourSideTransform(TPZGeoElSide &neighbour);
	
	/** @brief Compute the transformation between the master element space of one side
	 of an element to the master element space of a higher dimension side
	 */
	TPZTransform SideToSideTransform(TPZGeoElSide &higherdimensionside);
	
	/** Returns the number of fathers that can be followed*/
	//int Level();
	/** @brief Returns a pointer to the elementside referenced by the geometric elementside*/
	TPZCompElSide Reference() const;
	/** @brief Return 1 if the element has subelements along side*/
	int HasSubElement();
	
	/* @brief Return the number of nodes for a particular side*/
	int NSideNodes() const;
	
	/** @brief Returns the index of the nodenum node of side*/
	int SideNodeIndex(int nodenum) const;
	std::set<int> SideNodeIndexes();
	
	/** @brief Returns the index of the local nodenum node of side*/
	int SideNodeLocIndex(int nodenum) const;
	
	/** @brief Returns 1 if neighbour is a neighbour of the element along side*/
	int NeighbourExists(const TPZGeoElSide &neighbour) const;
	
	/** @brief Pushes all connected computational elements which have higher dimension than the current element/side
     
	 if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
	 */
	void HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated);
	
	
	/** Returns all connected computational elements which have same dimension and level higher to the current element
     
	 if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	//void HigherLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
	
	/** @brief Returns all connected computational elements which have level higher to the current element
     if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void HigherLevelCompElementList2(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
	/** @brief Returns all connected computational elements to the current element
     
	 if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void ConnectedCompElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);
	/** @brief Returns all connected computational elements which have equal level to the current element
     
	 This method will not put this on the stack
     if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,	int onlyinterpolated, int removeduplicates);
	//void Dim0EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);
	
	/** Returns all connected computational elements which have level lower to the current element
     if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
     
	 This method will only return element/sides of dimension > 0*/
	//TPZCompElSide LowerLevelCompElementList(int onlyinterpolated);
	/**Retorna o elem. computacional de nivel menor (elemento grande) ao qual estou restrito*/
	//TPZCompElSide TPZGeoElSide::Dim0LowerLevelCompElement(int onlyinterpolated);
	
	/** Jorge 13/01/2000 */
	/** Retorna todos os subelementos computacionais do atual elemento geometrico
	 * ao longo do lado fSide
	 */
	//void GetCompSubElements(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated,int removeduplicates);
	
};

/** @brief Overload operator << to print geometric element side data */
std::ostream  &operator << (std::ostream & out,const TPZGeoElSide &geoside);

#include "pzgeoel.h"

inline TPZGeoElSide TPZGeoElSide::Neighbour() const {
	return fGeoEl ? fGeoEl->Neighbour(fSide) : TPZGeoElSide();
}

inline void TPZGeoElSide::AllNeighbours(TPZStack<TPZGeoElSide> &allneigh) {
	TPZGeoElSide neigh = Neighbour();
#ifndef NODEBUG
	if(! Exists() || ! neigh.Exists()) 
    {
		std::cout << "TPZGeoElSide AllNeighbours inconsistent\n";
		return;
    }
#endif
	while(neigh != *this)
    {
		allneigh.Push(neigh);
		neigh = neigh.Neighbour();
    }
}


/** Implementing TPZGeoElSideIndex methods */

inline TPZGeoElSideIndex::~TPZGeoElSideIndex(){
	//nothing to be done here
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(){
    this->fSide = -1;
    this->fGeoElIndex = -1;
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(TPZGeoEl *gel,int side){  
    if (gel) this->fGeoElIndex = gel->Index();
    else this->fGeoElIndex = -1;
    this->fSide = side;
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(int gelindex,int side){  
    this->fGeoElIndex = gelindex;
    this->fSide = side;
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(const TPZGeoElSide &side){
    TPZGeoEl * gel = side.Element();
    if (gel) this->fGeoElIndex = gel->Index();
    else this->fGeoElIndex = -1;
    this->fSide = side.Side();
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(const TPZGeoElSideIndex &cp){
    this->operator =(cp);
}

inline TPZGeoElSideIndex * TPZGeoElSideIndex::Clone(){
    return new TPZGeoElSideIndex(*this);
}

inline TPZGeoElSideIndex & TPZGeoElSideIndex::operator= (const TPZGeoElSideIndex &A ){
    this->fGeoElIndex = A.fGeoElIndex;
    this->fSide = A.fSide;
    return *this;
}

inline int TPZGeoElSideIndex::Side() const{
	// #ifdef DEBUG
	//     if (this->fSide == -1 || this->fGeoElIndex == -1){
	//       PZError << "Error at " << __PRETTY_FUNCTION__ << " - TPZGeoElSideIndex not initialized is NULL\n";
	//     }
	// #endif
    if (this->fSide == -1 || this->fGeoElIndex == -1){
		return -1;
    }
    return this->fSide;
}

inline void TPZGeoElSideIndex::SetSide(int side){
    this->fSide = side;
}  

inline TPZGeoEl *TPZGeoElSideIndex::Element(const TPZGeoMesh *mesh) const{
    if (this->fSide == -1 || this->fGeoElIndex == -1){
		return NULL;
    }
    return mesh->ElementVec()[this->fGeoElIndex];
}

inline void TPZGeoElSideIndex::SetElement(TPZGeoEl* geoel){
    if (geoel) this->fGeoElIndex = geoel->Index();
    else this->fGeoElIndex = -1;
}

inline int TPZGeoElSideIndex::ElementIndex() const{
    return this->fGeoElIndex;
}

inline void TPZGeoElSideIndex::SetElementIndex(int i){
    this->fGeoElIndex = i;
}

inline void TPZGeoElSideIndex::Read(TPZStream &buf){
    int side, index;
    buf.Read(&side, 1);
    buf.Read(&index, 1);
    this->fSide = side;
    this->fGeoElIndex = index;
}

inline void TPZGeoElSideIndex::Write(TPZStream &buf){
    int side = this->fSide;
    int index = this->fGeoElIndex;
    buf.Write(&side, 1);
    buf.Write(&index, 1);
}

#endif
