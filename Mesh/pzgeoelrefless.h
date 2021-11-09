/**
 * @file
 * @brief Contains declaration of TPZGeoElRefLess class which implements the mapping between the master element and deformed element.
 */

#ifndef PZGEOELREFLESS_H
#define PZGEOELREFLESS_H
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgeom_utility.h"

class TPZGeoElSide;
class TPZCompMesh;
class TPZCompEl;
template<class T,int N>
class TPZStack;

/**
 * Implements a generic geometric element class without h-refinement. \n
 * Its data structure is the vector of node indices and element neighbours \n
 * Implements a class which subelement is the clone of the element (i.e. the same nodes, material, but a self pointer)
 */
/**
 * @brief Implements the mapping between the master element and deformed element. \ref geometry "Geometry"
 * @author Philippe Devloo
 * @ingroup geometry
 * @since Dez 12, 2003.
 */
template <class TGeo>
class TPZGeoElRefLess : public TPZGeoEl  {
	//  int fSubElement;
protected:
	TGeo fGeo;
	//  int fNodeIndexes[TGeo::NNodes];
	TPZGeoElSideIndex fNeighbours[TGeo::NSides];
public:

virtual int ClassId() const override;

	virtual ~TPZGeoElRefLess();
	TPZGeoElRefLess();
	
	/** @brief Copy constructor */
	TPZGeoElRefLess(const TPZGeoElRefLess &gel);
	
	virtual TPZGeoEl *Clone(TPZGeoMesh &dest) const override
	{
		return new TPZGeoElRefLess(dest,*this);
	}
	
	/** @brief Copy constructor with elements in different meshes */
	TPZGeoElRefLess(TPZGeoMesh &DestMesh, const TPZGeoElRefLess &cp);
    
    //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
	
	/**
	 * @brief Copy constructor with elements in different meshes. The clone mesh is
	 * a patch mesh.
	 *
	 * Therefore there are the requirement to map the nodes
	 * between these meshses.
	 * @param DestMesh destination patch mesh
	 * @param cp element to be copied
	 * @param gl2lcNdMap map of the node indexes between original and clone mesh
	 * @param gl2lcElMap map of the element indexes between original and clone mesh
	 */
	TPZGeoElRefLess(TPZGeoMesh &DestMesh,
					const TPZGeoElRefLess &cp,
					std::map<int64_t,int64_t> & gl2lcNdMap,
					std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZGeoEl *ClonePatchEl(TPZGeoMesh &destmesh,std::map<int64_t,int64_t> & gl2lcNdMap,
						   std::map<int64_t,int64_t> & gl2lcElMap) const  override
	{
		return new TPZGeoElRefLess(destmesh,*this,gl2lcNdMap, gl2lcElMap);
	}
	
	TPZGeoElRefLess(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
	
	TPZGeoElRefLess(TGeo &Geo, int matind,TPZGeoMesh &mesh);
	
	TPZGeoElRefLess(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh);
    
	TPZGeoElRefLess(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int64_t &index);
	
	virtual void Read(TPZStream &str, void *context) override;
	
	virtual void Write(TPZStream &str, int withclassid) const override;
	
	virtual void Initialize() override
	{
		fGeo.Initialize(this);
	}
	
	
	static int main_refless();
	
	/** @brief Divides the element and puts the resulting elements in the vector */
	virtual void Divide(TPZVec < TPZGeoEl * > & pv) override {
		DebugStop();
	}
	
	/** @brief Returns 1 if the element has subelements along side*/
	virtual  int HasSubElement() const  override {return 0;}//fSubEl[0]!=0;}
	
	/** @brief Returns a pointer to the neighbour and the neighbourside along side of the current element */
	virtual  TPZGeoElSide Neighbour(int side) override {
#ifdef PZDEBUG
        if (fNeighbours[side] < 0 || fNeighbours[side] >= this->Mesh()->NElements()) {
            DebugStop();
        }
#endif
        return TPZGeoElSide(fNeighbours[side],this->Mesh());
    }
	
	/** @brief Returns the neighbour index for a given side*/
	virtual int64_t NeighbourIndex(int side) const override{
#ifdef PZDEBUG
        	if (fNeighbours[side] < 0 || fNeighbours[side] >= this->Mesh()->NElements()) {
            		DebugStop();
        	}
#endif
        	return this->fNeighbours[side].ElementIndex();
	}
	
	virtual  int64_t NodeIndex(int node) const override;
	
    inline void CornerCoordinates(TPZFMatrix<REAL> &coord) const;
	//HDiv
    
//    virtual void Directions(int side, TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions, TPZVec<int> &vectorsides) override;

    virtual void HDivDirectionsMaster(TPZFMatrix<REAL> &directions) override;
    
    virtual void HDivDirections(TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions) override;
    

    virtual void HDivDirections(TPZVec<REAL> &pt, TPZFMatrix<Fad<REAL> > &directions) override;

	//virtual void VecHdiv(TPZFMatrix<REAL> &normalvec ,TPZVec<int> &sidevector) override;

	
	/** @brief Compute the permutation for an HDiv side */
	virtual void HDivPermutation(int side, TPZVec<int> &permutegather) override;
	
	
	/** @brief Fill in the data structure for the neighbouring information*/
	virtual  void SetNeighbour(int side,const TPZGeoElSide &neighbour) override {
		fNeighbours[side]=neighbour;
	}
	
	virtual void Print(std::ostream &out) override
	{
		TPZGeoEl::Print(out);
		out << "fGeo:\n";
		fGeo.Print(out);
	}
	/** @brief Prints topological information of: TGeo (TPZGeoCube, TPZGeoPrism, TPZGeoQuad, ...) */
	virtual void PrintTopologicalInfo(std::ostream &out) override
	{
		out << "Geo Element - fId " << fId << "\t Type " << fGeo.TypeName() << "\t";
		fGeo.Print(out);
		int i;
		out << std::endl << "\t";
		for (i = 0;i < NNodes();i++) out << NodePtr(i)->Id() << " ";
	}
	
	virtual  int64_t SideNodeIndex(int side,int node) const override;
	
	virtual  int SideNodeLocIndex(int side,int node) const override;
	
	/** @brief Flags the side as defined, this means no neighbouring element was found */
	virtual  void SetSideDefined(int side)  override { fNeighbours[side] = TPZGeoElSide(this,side); }
	
	virtual  void SetSubElement(int id, TPZGeoEl *el) override;


    /**
     * This method gets the ith valid permutation of its topology
     * @param i number of the permutation to get
     * @param permutation vector contained the permuted sides
     */
    void GetPermutation(const int& i, TPZVec<int> &permutation) const override;
	/**
	 * @brief Creates an integration rule for the topology of the corresponding side
	 * and able to integrate a polynom of order exactly
	 */
	virtual  TPZIntPoints * CreateSideIntegrationRule(int side, int order) override;
	
	/** @brief Returns the type of the element acording to the definition in pzeltype.h */
	virtual  MElementType Type() const override {
		return TGeo::Type();
	}
	
	/** @brief Returns the type of the element acording to the definition in pzeltype.h */
	virtual  MElementType Type(int side) const override {
		return TGeo::Type(side);
	}
	
	/** @brief Returns the type of the element as a string */
    virtual std::string TypeName() const override
    {
		return fGeo.TypeName();
    }
	
	/** @brief Returns the number of nodes of the element*/
	virtual  int NNodes() const override;
	
	/** @brief Returns the number of corner nodes of the element*/
	virtual  int NCornerNodes() const override;
	
	/** @brief Returns the number of sides of the element*/
	virtual  int NSides() const override;

    /** @brief Returns the number of sides of the element of a given dimension */
    virtual  int NSides(int dim) const override;

    /** @brief Returns the first side of the element of a given dimension */
    virtual  int FirstSide(int dim) const override;

    /**
     * Get the number of valid permutations among the element nodes
     * @return
     */
    virtual int NPermutations() const override;
	
	/** @brief Returns the local node number of the node "node" along side "side" */
	virtual  int SideNodeLocId(int side, int node) const;
	
	/** @brief Volume of the master element*/
	virtual  REAL RefElVolume() override;
	
	/** @brief Returns the number of nodes for a particular side*/
	virtual  int NSideNodes(int side) const override;
	
	/** @brief Returns the midside node index along a side of the element*/
	virtual  void MidSideNodeIndex(int side,int64_t &index) const override;
	
	/** @brief Returns 1 if the side has not been defined by buildconnectivity */
	/** After construction the side is undefined. The buildconnectivity method
	 * loops over all elements and tries to identify neighbours along their uninitialized sides */
	virtual  int SideIsUndefined(int side) override;
	
	/** @brief Returns the number of subelements of the element independent of the fact hether the element has already been refined or not */
	virtual  int NSubElements() const override;
	
	/** @brief Returns the number of subelements of the same dimension of the element at the side*/
	virtual  int NSideSubElements(int side) const override;
	
	/**
	 * @brief Method which creates a computational boundary condition element based
	 * on the current geometric element, a side and a boundary condition number
	 */
	virtual  TPZGeoEl *CreateBCGeoEl(int side, int bc) override;
	
	/**
	 * @brief Method which creates a blend geometrical boundary condition element
	 * based on the current geometric element, a side and a boundary condition index
	 */
	virtual TPZGeoEl *CreateBCGeoBlendEl(int side, int bc);

    /**
     * @brief Method which creates a blend geometrical boundary condition element
     * based on the current geometric element, a side and a boundary condition index
     */
    virtual TPZGeoEl *CreateBCGeoBlendEl(int side, int bc, const TPZVec<int> &mapsides);

	/** @brief Creates a geometric element according to the type of the father element */
	virtual TPZGeoEl *CreateGeoElement(MElementType type,
									   TPZVec<int64_t>& nodeindexes,
									   int matid,
									   int64_t& index) override;
	
	/** @brief Initializes the node i of the element*/
	virtual  void SetNodeIndex(int i,int64_t nodeindex) override;
	
	/**
	 * @brief compute the transformation between the master element space of one side
	 * of an element to the master element space of a higher dimension side
	 */
	virtual  TPZTransform<> SideToSideTransform(int sidefrom,int sideto) override;
	
	/** @brief Returns a pointer to the subelement is*/
	virtual  TPZGeoEl *SubElement(int is) const override;
	
	/** @brief Return the dimension of side*/
	virtual  int SideDimension(int side) const override;
	
	/** @brief Returns the dimension of the element*/
	virtual int Dimension() const override;
	
	virtual  TPZGeoElSide HigherDimensionSides(int side,int targetdimension);
	
	virtual  void AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides) override;
	
	virtual  void LowerDimensionSides(int side,TPZStack<int> &smallsides) const override;
	
	/** @brief Accumulates the transformation of the jacobian which maps the current
     master element space into the space of the master element of the father*/
	virtual  void BuildTransform(int side, TPZGeoEl *father,TPZTransform<> &t);
	
	virtual  TPZTransform<> BuildTransform2(int side, TPZGeoEl *father,TPZTransform<> &t) override;
	
    /** @brief Returns the coordinate in real space of the point coordinate in the master element space*/
    virtual  void X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result) const override;
    
    /** @brief Return the gradient of the transformation at the point */
    virtual void GradX(TPZVec<REAL> &coordinate, TPZFMatrix<REAL> &gradx) const override;
    
    /** @brief Returns the coordinate in real space of the point coordinate in the master element space*/
    virtual  void X(TPZVec<Fad<REAL> > &coordinate,TPZVec<Fad<REAL> > &result) const override;
    
    /** @brief Return the gradient of the transformation at the point */
    virtual void GradX(TPZVec<Fad<REAL> > &coordinate, TPZFMatrix<Fad<REAL> > &gradx) const override;

	virtual bool IsLinearMapping( int side) const override;
	virtual bool IsGeoBlendEl() const override;

    /**
     * If the element is a TPZGeoBlend element, this method will ensure that if the side side is connected to the
     * element with index index, its blend connectivity is erased (for instance, this element may have been deleted).
     * If it is not a TPZGeoBlend element, the method will just return false.
     * @param side side in which to seek for connectivities
     * @param index index of the element that will be disconnected from this
     * @return true if the element is a TPZGeoBlend element.
     */
    bool ResetBlendConnectivity(const int64_t &side, const int64_t &index) override;
	TGeo &Geom() { return fGeo; }
	
	virtual  TPZTransform<> GetTransform(int side,int son) override;
	
	/** @brief It returns the coordinates of the center of the side of the element */
	virtual void CenterPoint(int side, TPZVec<REAL> &masscent) const override;
	
	virtual TPZGeoElSide Father2(int side) const override;
	
	virtual int FatherSide(int side, int son)  override {
		return side;
	}
	
	virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const override;
	
	virtual void ResetSubElements() override {
		DebugStop();
	}
	
	virtual void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<> &trans) override
	{
		Geom().SetNeighbourInfo(side,neigh,trans);
	}
	
    /** @brief Generates a random point in the master domain */
    virtual void RandomPoint(TPZVec<REAL> &pt) override
    {
        Geom().RandomPoint(pt);
    }
    
	/** @brief Verifies if the parametric point pt is in the element parametric domain */
	virtual bool IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol = 1e-6) override;
	
	/**
	 * @brief Ortogonal projection from given qsi to a qsiInDomain (all in the element parametric domain)
	 * @return Returns the side where the point was projected.
     * @note If ortogonal projection to any sides of element results in a qsi outside domain,
     *       this method will returns the nearest node side.
	 * @note Observe that if the point is already in the parametric domain, the method will return \f$ NSides() - 1 \f$
	 */
	virtual int ProjectInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain) override;
    
    /**
	 * @brief Projection from given qsi to a qsiInDomain (in the element boundary) using bissection method from given qsi to element center.
	 * @return Returns the side where the point was projected.
	 * @note Observe that if the point is already in the parametric domain, the method will return \f$ NSides() - 1 \f$
	 */
    virtual int ProjectBissectionInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain) override;
};

template<class TGeo>
inline
bool TPZGeoElRefLess<TGeo>::IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol){
	const bool result = fGeo.IsInParametricDomain(pt,tol);
	return result;
}

template<class TGeo>
inline
int TPZGeoElRefLess<TGeo>::ProjectInParametricDomain(TPZVec<REAL> &pt, TPZVec<REAL> &ptInDomain){
    const int side = ::ProjectInParametricDomain<TGeo>(pt, ptInDomain);
	return side;
}

template<class TGeo>
inline
int TPZGeoElRefLess<TGeo>::ProjectBissectionInParametricDomain(TPZVec<REAL> &pt, TPZVec<REAL> &ptInDomain){
    const int side = ::ProjectBissectionInParametricDomain<TGeo>(pt, ptInDomain);
	return side;
}

template<class TGeo>
int TPZGeoElRefLess<TGeo>::ClassId() const{
    return Hash("TPZGeoElRefLess") ^ TPZGeoEl::ClassId() << 1 ^ TGeo().ClassId() << 2;
}

//template<class TGeo>
//inline
//void TPZGeoElRefLess<TGeo>::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    fGeo.ParametricDomainNodeCoord(node, nodeCoord);
//}

#include "pzgeoelrefless.h.h"

#endif
