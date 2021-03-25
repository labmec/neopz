/**
 * @file
 * @brief Contains the implementation of the TPZGeoElement methods.
 */

#include "TPZGeoElement.h"
#include <sstream>
#include "pzlog.h"

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh) :
TPZRegisterClassId(&TPZGeoElement::ClassId),TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh, int64_t &index) :
TPZRegisterClassId(&TPZGeoElement::ClassId),TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh,index) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TGeo &geo,int matind,TPZGeoMesh &mesh) :
TPZRegisterClassId(&TPZGeoElement::ClassId),TPZGeoElRefLess<TGeo>(geo,matind,mesh) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
TPZRegisterClassId(&TPZGeoElement::ClassId),TPZGeoElRefLess<TGeo>(id,nodeindexes,matind,mesh) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement() : TPZRegisterClassId(&TPZGeoElement::ClassId),TPZGeoElRefLess<TGeo>() {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
void TPZGeoElement<TGeo,TRef>::SetSubElement(int id, TPZGeoEl *el) {
	
	if (id<0 || id >(TRef::NSubEl - 1)){
		PZError << "TPZGeoElement::Trying do define subelement :"
	    << id << "Max Allowed = " << TRef::NSubEl - 1 << std::endl;
		return;
	}
    if (el)
    {
        fSubEl[id] = el->Index();
    }
    else
    {
        fSubEl[id] = -1;
    }
	return;
}

template<class TGeo, class TRef>
REAL TPZGeoElement<TGeo,TRef>::RefElVolume(){
	return TGeo::RefElVolume();
}

template<class TGeo, class TRef>
void TPZGeoElement<TGeo,TRef>::MidSideNodeIndex(int side,int64_t &index) const
{
	TRef::MidSideNodeIndex(this,side,index);
}


template<class TGeo, class TRef>
int TPZGeoElement<TGeo,TRef>::NSubElements() const 
{
	return TRef::NSubEl;
}

template<class TGeo, class TRef>
int TPZGeoElement<TGeo,TRef>::NSideSubElements(int side) const
{
	return TRef::NSideSubElements(side);
}

template<class TGeo, class TRef>
TPZGeoEl *TPZGeoElement<TGeo,TRef>::SubElement(int is) const
{
	if(is<0 || is>(TRef::NSubEl - 1)){
		std::cout << "TPZGeoElement::SubElement index error is= " << is << std::endl;
	}
	if(fSubEl[is] == -1) return 0;
	return this->Mesh()->ElementVec()[fSubEl[is]];
}
/*!
 \fn TPZGeoElement::ResetSubElements()
 */
template<class TGeo, class TRef>
void TPZGeoElement<TGeo, TRef>::ResetSubElements() {
    for (unsigned int i = 0; i < NSubElements(); ++i) {
        TPZGeoEl *gel = SubElement(i);
        if (gel) {
            gel->SetFatherIndex(-1);
        }
        fSubEl[i] = -1;
    }
}

template<class TGeo, class TRef>
TPZGeoElSide TPZGeoElement<TGeo,TRef>::SideSubElement(int side,int position) {
	TPZStack<TPZGeoElSide> subs;
	TRef::GetSubElements(this,side,subs);
	return subs[position];
}

template<class TGeo, class TRef>
TPZTransform<> TPZGeoElement<TGeo,TRef>::GetTransform(int side,int son) {
	return TRef::GetTransform(side,son);
}

template<class TGeo, class TRef>
void TPZGeoElement<TGeo,TRef>::Divide(TPZVec<TPZGeoEl *> &pv){
	
	TRef::Divide(this,pv);
}

template<class TGeo, class TRef>
void TPZGeoElement<TGeo,TRef>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const
{
	
	TRef::GetSubElements(this,side,subel);
}

template<class TGeo, class TRef>
void TPZGeoElement<TGeo,TRef>::Read(TPZStream &buf, void *context) {
	TPZGeoElRefLess<TGeo>::Read(buf,context);
	buf.Read(fSubEl,TRef::NSubEl);
}

template<class TGeo, class TRef>
void TPZGeoElement<TGeo,TRef>::Write(TPZStream &buf, int withclassid) const{
	TPZGeoElRefLess<TGeo>::Write(buf,withclassid);
	buf.Write(fSubEl,TRef::NSubEl);
}

template<class TGeo, class TRef>
TPZGeoEl * TPZGeoElement<TGeo,TRef>::Clone(TPZGeoMesh &DestMesh) const{
    if(&DestMesh == this->Mesh())
    {
        return new TPZGeoElement<TGeo,TRef>(*this);
    }
    else
    {
        return new TPZGeoElement<TGeo,TRef>(DestMesh, *this);
    }
}//Clone method

template<class TGeo, class TRef>
TPZGeoEl * TPZGeoElement<TGeo,TRef>::ClonePatchEl(TPZGeoMesh &DestMesh,
												  std::map<int64_t,int64_t> & gl2lcNdMap,
												  std::map<int64_t,int64_t> & gl2lcElMap) const{
	return new TPZGeoElement<TGeo,TRef>(DestMesh, *this, gl2lcNdMap, gl2lcElMap);
}//Clone method


template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZGeoMesh &DestMesh, const TPZGeoElement &cp):
TPZRegisterClassId(&TPZGeoElement::ClassId),
TPZGeoElRefLess<TGeo>(DestMesh, cp){
	int i, n = TRef::NSubEl;
	for(i = 0; i < n; i++){
		this->fSubEl[i] = cp.fSubEl[i];
	}
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZGeoMesh &DestMesh,
										const TPZGeoElement &cp,
										std::map<int64_t,int64_t> &gl2lcNdMap,
										std::map<int64_t,int64_t> &gl2lcElMap):
TPZRegisterClassId(&TPZGeoElement::ClassId),
TPZGeoElRefLess<TGeo>(DestMesh, cp, gl2lcNdMap, gl2lcElMap)
{
	int i, n = TRef::NSubEl;
	for(i = 0; i < n; i++)
	{
		if (cp.fSubEl[i] == -1)
		{
			this->fSubEl[i]=-1;
			continue;
		}
		if (gl2lcElMap.find(cp.fSubEl[i]) == gl2lcElMap.end())
		{
#ifdef PZ_LOG
            TPZLogger logger("pz.mesh.tpzgeoelement");

            if (logger.isDebugEnabled()) {
              std::stringstream sout;
              sout << "ERROR in - " << __PRETTY_FUNCTION__
                   << " subelement index is not in map from original to "
                      "clone indexes! Son index = "
                   << fSubEl[i];
              LOGPZ_ERROR(logger, sout.str().c_str());
            }
#endif
            DebugStop();
        }
		this->fSubEl[i] = gl2lcElMap[cp.fSubEl[i]];
	}
}
