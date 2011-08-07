/**
 * @file
 * @brief Contains the implementation of the TPZGeoElement methods.
 */
//
// C++ Interface: TPZGeoElement
//
// Description:
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "TPZGeoElement.h"
#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeoelement"));
#endif

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) :
TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) :
TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh,index) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TGeo &geo,int matind,TPZGeoMesh &mesh) :
TPZGeoElRefLess<TGeo>(geo,matind,mesh) {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
TPZGeoElRefLess<TGeo>(id,nodeindexes,matind,mesh) {
	
	//  int i,nnod = nodeindexes.NElements();
	//  if(nnod!=TGeo::NNodes) {
	//    PZError << "TPZGeoElement<TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
	//    return;
	//  }
	//
	//  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindexes[i];
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

/*
 template< class TGeo, class TRef >
 void TPZGeoElement< TGeo, TRef >::Initialize(TPZVec<int> &nodeindices ) {
 
 TPZGeoElRefLess<TGeo>::Initialize(nodeindices);
 for( int i = 0; i < TRef::NSubEl; i++ ){
 fSubEl[ i ] = 0;
 }
 }
 */

template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement() : TPZGeoElRefLess<TGeo>() {
	int i;
	for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}


template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::SetSubElement(int id, TPZGeoEl *el){
	
	if (id<0 || id >(TRef::NSubEl - 1)){
		PZError << "TPZGeoElement::Trying do define subelement :"
	    << id << "Max Allowed = " << TRef::NSubEl - 1 << std::endl;
		return;
	}
	fSubEl[id] = el->Index();
	return;
}


template<class TGeo, class TRef>
REAL
TPZGeoElement<TGeo,TRef>::RefElVolume(){
	return TGeo::RefElVolume();
}


template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::MidSideNodeIndex(int side,int &index){
	TRef::MidSideNodeIndex(this,side,index);
}


template<class TGeo, class TRef>
int
TPZGeoElement<TGeo,TRef>::NSubElements(){
	return TRef::NSubEl;
}


template<class TGeo, class TRef>
int
TPZGeoElement<TGeo,TRef>::NSideSubElements2(int side){
	return TRef::NSideSubElements(side);
}


template<class TGeo, class TRef>
TPZGeoEl *
TPZGeoElement<TGeo,TRef>::SubElement(int is){
	if(is<0 || is>(TRef::NSubEl - 1)){
		std::cout << "TPZGeoElement::SubElement index error is= " << is << std::endl;;
	}
	if(fSubEl[is] == -1) return 0;
	return this->Mesh()->ElementVec()[fSubEl[is]];
}


/*!
 \fn TPZGeoElement::ResetSubElements()
 */
template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::ResetSubElements()
{
	int is;
	for (is=0;is<NSubElements();is++)
	{
		fSubEl[is] = -1;
	}
}


template<class TGeo, class TRef>
TPZGeoElSide
TPZGeoElement<TGeo,TRef>::SideSubElement(int side,int position){
	TPZStack<TPZGeoElSide> subs;
	TRef::GetSubElements(this,side,subs);
	return subs[position];
}


template<class TGeo, class TRef>
TPZTransform
TPZGeoElement<TGeo,TRef>::GetTransform(int side,int son){
	
	return TRef::GetTransform(side,son);
}




template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::Divide(TPZVec<TPZGeoEl *> &pv){
	
	TRef::Divide(this,pv);
}

template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
	
	TRef::GetSubElements(this,side,subel);
}

template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::Read(TPZStream &buf, void *context) {
	TPZGeoElRefLess<TGeo>::Read(buf,context);
	buf.Read(fSubEl,TRef::NSubEl);
}

template<class TGeo, class TRef>
void
TPZGeoElement<TGeo,TRef>::Write(TPZStream &buf, int withclassid) {
	TPZGeoElRefLess<TGeo>::Write(buf,withclassid);
	buf.Write(fSubEl,TRef::NSubEl);
}

template<class TGeo, class TRef>
TPZGeoEl * TPZGeoElement<TGeo,TRef>::Clone(TPZGeoMesh &DestMesh) const{
	return new TPZGeoElement<TGeo,TRef>(DestMesh, *this);
}//Clone method


template<class TGeo, class TRef>
TPZGeoEl * TPZGeoElement<TGeo,TRef>::ClonePatchEl(TPZGeoMesh &DestMesh,
												  std::map<int,int> & gl2lcNdMap,
												  std::map<int,int> & gl2lcElMap) const{
	return new TPZGeoElement<TGeo,TRef>(DestMesh, *this, gl2lcNdMap, gl2lcElMap);
}//Clone method


template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZGeoMesh &DestMesh, const TPZGeoElement &cp):
TPZGeoElRefLess<TGeo>(DestMesh, cp){
	int i, n = TRef::NSubEl;
	for(i = 0; i < n; i++){
		this->fSubEl[i] = cp.fSubEl[i];
	}
}


template<class TGeo, class TRef>
TPZGeoElement<TGeo,TRef>::TPZGeoElement(TPZGeoMesh &DestMesh,
										const TPZGeoElement &cp,
										std::map<int,int> &gl2lcNdMap,
										std::map<int,int> &gl2lcElMap):
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
			std::stringstream sout;
			sout << "ERROR in - " << __PRETTY_FUNCTION__
			<< " subelement index is not in map from original to clone indexes! Son index = "
			<<  fSubEl[i];
			LOGPZ_ERROR(logger,sout.str().c_str());
			exit(-1);
		}
		this->fSubEl[i] = gl2lcElMap[cp.fSubEl[i]];
	}
}

