//#include "pzelgpoint.h"
//#include "pzelg1d.h"
//#include "pzelgt2d.h"
//#include "pzelgq2d.h"
//#include "pzelgt3d.h"
//#include "pzelgpi3d.h"
//#include "pzelgpr3d.h"
//#include "pzelgc3d.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoElement.h"
//#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstream.h"
#include "pzmeshid.h"
//#include "pzstack.h"

//template<class TShape, class TGeo, class TRef>
//int TPZGeoElement<TShape,TGeo,TRef>::fTest;

template<class TShape, class TGeo, class TRef>
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) : 
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh) {
  
//  int i,nnod = nodeindices.NElements();
//  if(nnod!=TGeo::NNodes) {
//    PZError << "TPZGeoElement<TShape,TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
//    return;
//  }
  
//  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  int i;
  for(i=0;i<TGeo::NNodes;i++) fSubEl[i] = -1;
}

template<class TShape, class TGeo, class TRef>
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) : 
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh,index) {
  
//  int i,nnod = nodeindices.NElements();
//  if(nnod!=TGeo::NNodes) {
 //   PZError << "TPZGeoElement<TShape,TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
  //  return;
 // }
  
//  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  int i;
  for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template<class TShape, class TGeo, class TRef>
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoElRefLess<TShape,TGeo>(id,nodeindexes,matind,mesh) {

//  int i,nnod = nodeindexes.NElements();
//  if(nnod!=TGeo::NNodes) {
//    PZError << "TPZGeoElement<TShape,TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
//    return;
//  }
//
//  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindexes[i];
  int i;
  for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = -1;
}

template< class TShape, class TGeo, class TRef >
void TPZGeoElement< TShape, TGeo, TRef >::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index ) {
  
TPZGeoElRefLess<TShape,TGeo>::Initialize(nodeindices,matind,mesh,index);
  for( int i = 0; i < TRef::NSubEl; i++ ){
     fSubEl[ i ] = 0;
  }
}

template<class TShape, class TGeo, class TRef>
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement() : TPZGeoElRefLess<TShape,TGeo>() {
  int i;
  for(i=0;i<TGeo::NNodes;i++) fSubEl[i] = -1;
}


template<class TShape, class TGeo, class TRef>
void 
TPZGeoElement<TShape,TGeo,TRef>::SetSubElement(int id, TPZGeoEl *el){

  if (id<0 || id >(TRef::NSubEl - 1)){
    PZError << "TPZGeoElement::Trying do define subelement :" 
	    << id << "Max Allowed = " << TRef::NSubEl - 1 << endl;
    return;
  }
  fSubEl[id] = el->Index();
  return;
}


template<class TShape, class TGeo, class TRef>
REAL
TPZGeoElement<TShape,TGeo,TRef>::RefElVolume(){
	return TShape::RefElVolume();
}


template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::MidSideNodeIndex(int side,int &index){
	TRef::MidSideNodeIndex(this,side,index);
}


template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSubElements(){
	return TRef::NSubEl;
}


template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSideSubElements2(int side){
	return TRef::NSideSubElements(side);
}


template<class TShape, class TGeo, class TRef>
TPZGeoEl *
TPZGeoElement<TShape,TGeo,TRef>::SubElement(int is){
  if(is<0 || is>(TRef::NSubEl - 1)){
    cout << "TPZGeoElement::SubElement index error is= " << is << endl;;
  }
  if(fSubEl[is] == -1) return 0;
  return this->Mesh()->ElementVec()[fSubEl[is]];
}

template<class TShape, class TGeo, class TRef>
TPZGeoElSide
TPZGeoElement<TShape,TGeo,TRef>::SideSubElement(int side,int position){
  TPZStack<TPZGeoElSide> subs;
  TRef::GetSubElements(this,side,subs);
  return subs[position];
}


template<class TShape, class TGeo, class TRef>
TPZTransform
TPZGeoElement<TShape,TGeo,TRef>::GetTransform(int side,int son){

	return TRef::GetTransform(side,son);
}




template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::Divide(TPZVec<TPZGeoEl *> &pv){

  TRef::Divide(this,pv);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

  TRef::GetSubElements(this,side,subel);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::Read(TPZStream &buf, void *context) {
  TPZGeoElRefLess<TShape,TGeo>::Read(buf,context);
  buf.Read(fSubEl,TRef::NSubEl);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::Write(TPZStream &buf, int withclassid) {
  TPZGeoElRefLess<TShape,TGeo>::Write(buf,withclassid);
  buf.Write(fSubEl,TRef::NSubEl);
}

/*
template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::ClassId() const {
  return -1;
}
*/
template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
template class TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>;
template class TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>;
template class TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>;
template class TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>;
template class TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>;
template class TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>;
template class TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>;


int
TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>::ClassId() const {
  return TPZFGEOELEMENTPOINTID;
}
int
TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::ClassId() const {
  return TPZFGEOELEMENTLINEARID;
}
int
TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::ClassId() const {
  return TPZFGEOELEMENTQUADID;
}
int
TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::ClassId() const {
  return TPZFGEOELEMENTRIANGLEID;
}
int
TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::ClassId() const {
  return TPZFGEOELEMENTCUBEID;
}
int
TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::ClassId() const {
  return TPZFGEOELEMENTPRISMID;
}
int
TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::ClassId() const {
  return TPZFGEOELEMENTTETRAID;
}
int
TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::ClassId() const {
  return TPZFGEOELEMENTPYRAMID;
}
