/**
 * @file
 * @brief Contains declaration of TPZGeoElement class which implements a generic geometric element with a uniform refinement pattern.
 */

#ifndef TPZGEOELEMENTH
#define TPZGEOELEMENTH

#include "pzgeoelrefless.h"

class TPZGeoElSide;
class TPZCompMesh;
class TPZCompEl;
template<class T,int N>
class TPZStack;

/**
 * @ingroup geometry
 * @brief Implements a generic geometric element with a uniform refinement pattern. \ref geometry "Geometry"
 */
template <class TGeo, class TRef>
class TPZGeoElement : public TPZGeoElRefLess<TGeo> {
	
	long fSubEl[TRef::NSubEl];
public:
	typedef TGeo TGeoLoc;
	
public:
	/** @brief Default constructor */
	TPZGeoElement();
	/** @brief Constructor from node indexes and id given */
	TPZGeoElement(long id,TPZVec<long> &nodeindexes,int matind,TPZGeoMesh &mesh);
	/** @brief Constructor from node indexes */
	TPZGeoElement(TPZVec<long> &nodeindices,int matind,TPZGeoMesh &mesh);
	/** @brief Constructor with topology given */
	TPZGeoElement(TGeo &geo, int matind, TPZGeoMesh &mesh);
	/** @brief Constructor from node indexes and return the index of the new geometric element */
	TPZGeoElement(TPZVec<long> &nodeindices,int matind,TPZGeoMesh &mesh,long &index);
	
	/** @brief Copy constructor */
	TPZGeoElement(TPZGeoMesh &DestMesh, const TPZGeoElement &cp);
	
	/**
	 * @brief Clone constructor for patch meshes
	 * @param DestMesh destination patch mesh
	 * @param cp element to be copied
	 * @param gl2lcNdIdx map between the original node indexes and patch node indexes
	 * @param gl2lcElIdx map between the original element indexes and patch element index
	 */
	TPZGeoElement ( TPZGeoMesh &DestMesh,
				   const TPZGeoElement &cp,
				   std::map<long,long> &gl2lcNdIdx,
				   std::map<long,long> &gl2lcElIdx );
	/** @brief Default destructor */
	virtual ~TPZGeoElement(){};
	
	virtual int ClassId() const;
	
	virtual void Read(TPZStream &str, void *context);
	
	virtual void Write(TPZStream &str, int withclassid);
	
	virtual TPZGeoEl * Clone(TPZGeoMesh &DestMesh) const;
	
	/**
	 * @see class TPZGeoEl
	 */
	virtual TPZGeoEl * ClonePatchEl(TPZGeoMesh &DestMesh,
									std::map<long,long> &gl2lcNdIdx,
									std::map<long,long> &gl2lcElIdx) const;
	
	/** @brief Returns 1 if the element has subelements. */
	int HasSubElement() const {return fSubEl[0]!=-1;}
	
	void SetSubElement(int id, TPZGeoEl *el);
	
	
	/** @brief Volume of the master element*/
	REAL RefElVolume();
	
	/** @brief Returns the midside node index along a side of the element*/
	void MidSideNodeIndex(int side,long &index) const;
	
	/** @brief Returns the number of subelements of the element independent of the fact hether the element has already been refined or not */
	int NSubElements() const;
	
	/** @brief Returns the number of subelements as returned by GetSubElements2(side) */
	int NSideSubElements(int side) const;
	
	/** @brief Returns a pointer to the subelement is*/
	TPZGeoEl *SubElement(int is) const;
	
	/** @brief Return a pointer and a side of the subelement of the element at the side
     and the indicated position. position = 0 indicate first subelement, ...*/
	TPZGeoElSide SideSubElement(int side,int position);
	
	TPZTransform GetTransform(int side,int son);
	
	virtual int FatherSide(int side, int son) {
		return TRef::FatherSide(side,son);
	}
	
	/** @brief Divides the element and puts the resulting elements in the vector*/
	virtual void Divide(TPZVec<TPZGeoEl *> &pv);
	
	virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const;
	
	virtual void ResetSubElements();
	
};

#endif
