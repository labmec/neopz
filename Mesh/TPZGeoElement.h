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
	
	int64_t fSubEl[TRef::NSubEl];
public:
	typedef TGeo TGeoLoc;
	
public:
	/** @brief Default constructor */
	TPZGeoElement();
	/** @brief Constructor from node indexes and id given */
	TPZGeoElement(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
	/** @brief Constructor from node indexes */
	TPZGeoElement(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh);
	/** @brief Constructor with topology given */
	TPZGeoElement(TGeo &geo, int matind, TPZGeoMesh &mesh);
	/** @brief Constructor from node indexes and return the index of the new geometric element */
	TPZGeoElement(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int64_t &index);
	
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
				   std::map<int64_t,int64_t> &gl2lcNdIdx,
				   std::map<int64_t,int64_t> &gl2lcElIdx );
	/** @brief Default destructor */
	virtual ~TPZGeoElement(){};
	
	public:
int ClassId() const override;

	
	void Read(TPZStream &str, void *context) override;
	
	void Write(TPZStream &str, int withclassid) const override;
	
	virtual TPZGeoEl * Clone(TPZGeoMesh &DestMesh) const override;
	
	/**
	 * @see class TPZGeoEl
	 */
	virtual TPZGeoEl * ClonePatchEl(TPZGeoMesh &DestMesh,
									std::map<int64_t,int64_t> &gl2lcNdIdx,
									std::map<int64_t,int64_t> &gl2lcElIdx) const override;
	
	/** @brief Returns 1 if the element has subelements. */
	int HasSubElement() const override {return fSubEl[0]!=-1;}
	
	void SetSubElement(int id, TPZGeoEl *el) override;
	
	
	/** @brief Volume of the master element*/
	REAL RefElVolume() override;
	
	/** @brief Returns the midside node index along a side of the element*/
	void MidSideNodeIndex(int side,int64_t &index) const override;
	
	/** @brief Returns the number of subelements of the element independent of the fact hether the element has already been refined or not */
	int NSubElements() const override;
	
	/** @brief Returns the number of subelements as returned by GetSubElements2(side) */
	int NSideSubElements(int side) const override;
	
	/** @brief Returns a pointer to the subelement is*/
	TPZGeoEl *SubElement(int is) const override;
	
	/** @brief Return a pointer and a side of the subelement of the element at the side
     and the indicated position. position = 0 indicate first subelement, ...*/
	TPZGeoElSide SideSubElement(int side,int position);
	
	TPZTransform<> GetTransform(int side,int son) override;
	
	virtual int FatherSide(int side, int son)  override {
		return TRef::FatherSide(side,son);
	}
	
	/** @brief Divides the element and puts the resulting elements in the vector*/
	virtual void Divide(TPZVec<TPZGeoEl *> &pv) override;
	
	virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const override;
	
	virtual void ResetSubElements() override;
    
    /** @brief Creates a geometric element according to the type of the father element */
    virtual TPZGeoEl *CreateGeoElement(MElementType type,
                                       TPZVec<int64_t>& nodeindexes,
                                       int matid,
                                       int64_t& index) override
    {
        return this->Mesh()->CreateGeoElement(type,nodeindexes,matid,index,0);
    }
    

	
};

template <class TGeo, class TRef>
int TPZGeoElement<TGeo, TRef> ::ClassId() const{
    return Hash("TPZGeoElement") ^ TPZGeoElRefLess<TGeo>::ClassId() << 1 ^ TRef().ClassId() << 2;
}

#include "TPZGeoElement.h.h"
#endif
