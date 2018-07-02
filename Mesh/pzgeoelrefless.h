/**
 * @file
 * @brief Contains declaration of TPZGeoElRefLess class which implements the mapping between the master element and deformed element.
 */

#ifndef PZGEOELREFLESS_H
#define PZGEOELREFLESS_H

#include "pzgeoel.h"
#include "pzgeoelside.h"

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
virtual int ClassId() const;

	virtual ~TPZGeoElRefLess();
	TPZGeoElRefLess();
	
	/** @brief Copy constructor */
	TPZGeoElRefLess(const TPZGeoElRefLess &gel);
	
	virtual TPZGeoEl *Clone(TPZGeoMesh &dest) const
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
						   std::map<int64_t,int64_t> & gl2lcElMap) const
	{
		return new TPZGeoElRefLess(destmesh,*this,gl2lcNdMap, gl2lcElMap);
	}
	
	TPZGeoElRefLess(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
	
	TPZGeoElRefLess(TGeo &Geo, int matind,TPZGeoMesh &mesh);
	
	TPZGeoElRefLess(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh);
    
	TPZGeoElRefLess(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int64_t &index);
	
	virtual void Read(TPZStream &str, void *context);
	
	virtual void Write(TPZStream &str, int withclassid) const;
	
	virtual void Initialize()
	{
		fGeo.Initialize(this);
	}
	
	
	static int main_refless();
	
	/** @brief Divides the element and puts the resulting elements in the vector */
	virtual void Divide(TPZVec < TPZGeoEl * > & pv){
		DebugStop();
	}
	
	/** @brief Returns 1 if the element has subelements along side*/
	virtual  int HasSubElement() const {return 0;}//fSubEl[0]!=0;}
	
	/** @brief Returns a pointer to the neighbour and the neighbourside along side of the current element */
	virtual  TPZGeoElSide Neighbour(int side) {
#ifdef PZDEBUG
        if (fNeighbours[side] < 0 || fNeighbours[side] >= this->Mesh()->NElements()) {
            DebugStop();
        }
#endif
        return TPZGeoElSide(fNeighbours[side],this->Mesh());
    }
	
	virtual  int64_t NodeIndex(int node) const;
	
	//HDiv
    
    virtual void Directions(int side, TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions, TPZVec<int> &vectorsides);
    
    virtual void Directions(TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions, int ConstrainedFace = -1);
    
	virtual void VecHdiv(TPZFMatrix<REAL> &normalvec ,TPZVec<int> &sidevector);
	
	/** @brief Compute the permutation for an HDiv side */
	virtual void HDivPermutation(int side, TPZVec<int> &permutegather);
	
	
	/** @brief Fill in the data structure for the neighbouring information*/
	virtual  void SetNeighbour(int side,const TPZGeoElSide &neighbour){
		fNeighbours[side]=neighbour;
	}
	
	virtual void Print(std::ostream &out)
	{
		TPZGeoEl::Print(out);
		out << "fGeo:\n";
		fGeo.Print(out);
	}
	/** @brief Prints topological information of: TGeo (TPZGeoCube, TPZGeoPrism, TPZGeoQuad, ...) */
	virtual void PrintTopologicalInfo(std::ostream &out)
	{
		out << "Geo Element - fId " << fId << "\t Type " << fGeo.TypeName() << "\t";
		fGeo.Print(out);
		int i;
		out << std::endl << "\t";
		for (i = 0;i < NNodes();i++) out << NodePtr(i)->Id() << " ";
	}
	
	virtual  int64_t SideNodeIndex(int side,int node) const;
	
	virtual  int SideNodeLocIndex(int side,int node) const;
	
	/** @brief Flags the side as defined, this means no neighbouring element was found */
	virtual  void SetSideDefined(int side) { fNeighbours[side] = TPZGeoElSide(this,side); }
	
	virtual  void SetSubElement(int id, TPZGeoEl *el);
	
	/**
	 * @brief Creates an integration rule for the topology of the corresponding side
	 * and able to integrate a polynom of order exactly
	 */
	virtual  TPZIntPoints * CreateSideIntegrationRule(int side, int order);
	
	/** @brief Returns the type of the element acording to the definition in pzeltype.h */
	virtual  MElementType Type() const {
		return TGeo::Type();
	}
	
	/** @brief Returns the type of the element acording to the definition in pzeltype.h */
	virtual  MElementType Type(int side) const {
		return TGeo::Type(side);
	}
	
	/** @brief Returns the type of the element as a string */
    virtual std::string TypeName() const
    {
		return fGeo.TypeName();
    }
	
	/** @brief Returns the number of nodes of the element*/
	virtual  int NNodes() const;
	
	/** @brief Returns the number of corner nodes of the element*/
	virtual  int NCornerNodes() const;
	
	/** @brief Returns the number of connectivities of the element*/
	virtual  int NSides() const;
	
	/** @brief Returns the local node number of the node "node" along side "side" */
	virtual  int SideNodeLocId(int side, int node) const;
	
	/** @brief Volume of the master element*/
	virtual  REAL RefElVolume();
	
	/** @brief Returns the number of nodes for a particular side*/
	virtual  int NSideNodes(int side) const;
	
	/** @brief Returns the midside node index along a side of the element*/
	virtual  void MidSideNodeIndex(int side,int64_t &index) const;
	
	/** @brief Returns 1 if the side has not been defined by buildconnectivity */
	/** After construction the side is undefined. The buildconnectivity method
	 * loops over all elements and tries to identify neighbours along their uninitialized sides */
	virtual  int SideIsUndefined(int side);
	
	/** @brief Returns the number of subelements of the element independent of the fact hether the element has already been refined or not */
	virtual  int NSubElements() const;
	
	/** @brief Returns the number of subelements of the same dimension of the element at the side*/
	virtual  int NSideSubElements(int side) const;
	
	/**
	 * @brief Method which creates a computational boundary condition element based
	 * on the current geometric element, a side and a boundary condition number
	 */
	virtual  TPZGeoEl *CreateBCGeoEl(int side, int bc);
	
	/** @brief Creates a geometric element according to the type of the father element */
	virtual TPZGeoEl *CreateGeoElement(MElementType type,
									   TPZVec<int64_t>& nodeindexes,
									   int matid,
									   int64_t& index);
	
	/** @brief Initializes the node i of the element*/
	virtual  void SetNodeIndex(int i,int64_t nodeindex);
	
	/**
	 * @brief compute the transformation between the master element space of one side
	 * of an element to the master element space of a higher dimension side
	 */
	virtual  TPZTransform<> SideToSideTransform(int sidefrom,int sideto);
	
	/** @brief Returns a pointer to the subelement is*/
	virtual  TPZGeoEl *SubElement(int is) const;
	
	/** @brief Return the dimension of side*/
	virtual  int SideDimension(int side) const;
	
	/** @brief Returns the dimension of the element*/
	virtual int Dimension() const;
	
	virtual  TPZGeoElSide HigherDimensionSides(int side,int targetdimension);
	
	virtual  void AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides);
	
	virtual  void LowerDimensionSides(int side,TPZStack<int> &smallsides) const;
	
	/** @brief Accumulates the transformation of the jacobian which maps the current
     master element space into the space of the master element of the father*/
	virtual  void BuildTransform(int side, TPZGeoEl *father,TPZTransform<> &t);
	
	virtual  TPZTransform<> BuildTransform2(int side, TPZGeoEl *father,TPZTransform<> &t);
	
    /** @brief Returns the coordinate in real space of the point coordinate in the master element space*/
    virtual  void X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result) const;
    
    /** @brief Return the gradient of the transformation at the point */
    virtual void GradX(TPZVec<REAL> &coordinate, TPZFMatrix<REAL> &gradx) const ;
    
#ifdef _AUTODIFF
    /** @brief Returns the coordinate in real space of the point coordinate in the master element space*/
    virtual  void X(TPZVec<Fad<REAL> > &coordinate,TPZVec<Fad<REAL> > &result) const;
    
    /** @brief Return the gradient of the transformation at the point */
    virtual void GradX(TPZVec<Fad<REAL> > &coordinate, TPZFMatrix<Fad<REAL> > &gradx) const ;
#endif
    
	virtual bool IsLinearMapping( int side) const;
	virtual bool IsGeoBlendEl() const;
	TGeo &Geom() { return fGeo; }
	
	virtual  TPZTransform<> GetTransform(int side,int son);
	
	/** @brief It returns the coordinates of the center of the side of the element */
	virtual void CenterPoint(int side, TPZVec<REAL> &masscent) const;
	
	virtual TPZGeoElSide Father2(int side) const;
	
	virtual int FatherSide(int side, int son) {
		return side;
	}
	
	virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const;
	
	virtual void ResetSubElements(){
		DebugStop();
	}
	
	virtual void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<> &trans)
	{
		Geom().SetNeighbourInfo(side,neigh,trans);
	}
	
    /** @brief Generates a random point in the master domain */
    virtual void RandomPoint(TPZVec<REAL> &pt)
    {
        Geom().RandomPoint(pt);
    }
    
	/** @brief Verifies if the parametric point pt is in the element parametric domain */
	virtual bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);
	
	/**
	 * @brief Ortogonal projection from given qsi to a qsiInDomain (all in the element parametric domain)
	 * @return Returns the side where the point was projected.
     * @note If ortogonal projection to any sides of element results in a qsi outside domain,
     *       this method will returns the nearest node side.
	 * @note Observe that if the point is already in the parametric domain, the method will return \f$ NSides() - 1 \f$
	 */
	virtual int ProjectInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain);
    
    /**
	 * @brief Projection from given qsi to a qsiInDomain (in the element boundary) using bissection method from given qsi to element center.
	 * @return Returns the side where the point was projected.
	 * @note Observe that if the point is already in the parametric domain, the method will return \f$ NSides() - 1 \f$
	 */
    virtual int ProjectBissectionInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain);
};

template<class TGeo>
inline
bool TPZGeoElRefLess<TGeo>::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
	const bool result = fGeo.IsInParametricDomain(pt,tol);
	return result;
}

template<class TGeo>
inline
int TPZGeoElRefLess<TGeo>::ProjectInParametricDomain(TPZVec<REAL> &pt, TPZVec<REAL> &ptInDomain){
	const int side = fGeo.ProjectInParametricDomain(pt, ptInDomain);
	return side;
}

template<class TGeo>
inline
int TPZGeoElRefLess<TGeo>::ProjectBissectionInParametricDomain(TPZVec<REAL> &pt, TPZVec<REAL> &ptInDomain){
	const int side = fGeo.ProjectBissectionInParametricDomain(pt, ptInDomain);
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
