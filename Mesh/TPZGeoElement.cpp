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
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement() {
  int no;
  cout << "TPZGeoElement::TPZGeoElement constructor not implemented yet -> ";
  cin >> no;
}


template<class TShape, class TGeo, class TRef>
void 
TPZGeoElement<TShape,TGeo,TRef>::SetSubElement(int id, TPZGeoEl *el){

  if (id<0 || id >(TRef::NSubEl - 1)){
    PZError << "TPZGeoElement::Trying do define subelement :" 
	    << id << endl;
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
  return Mesh()->ElementVec()[fSubEl[is]];
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

template class TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>;
template class TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>;
template class TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>;
template class TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>;
template class TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>;
template class TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>;
template class TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>;
template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
/*
int newteste(){
TPZGeoElement <TPZShapeCube,TPZGeoCube,TPZRefCube> el1;
TPZGeoElement <TPZShapeLinear,TPZGeoLinear,TPZRefLinear> el2;
TPZGeoElement <TPZShapeQuad,TPZGeoQuad,TPZRefQuad> el3;
TPZGeoElement <TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle> el4;
TPZGeoElement <TPZShapePrism,TPZGeoPrism,TPZRefPrism> el5;
TPZGeoElement <TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra> el6;
TPZGeoElement <TPZShapePiram,TPZGeoPyramid,TPZRefPyramid> el7;
TPZGeoElement <TPZShapePoint,TPZGeoPoint,TPZRefPoint> el8;

return 0;
}
  
#include "pzelcpoint.h"
#include "pzelcq2d.h"
#include "pzelc1d.h"
#include "pzelct2d.h"
#include "pzelcc3d.h"
#include "pzelct3d.h"
#include "pzelcpr3d.h"
#include "pzelcpi3d.h"
#include "pzelctemp.h"

static TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoPoint,TPZShapePoint>(mesh,gel,index);
  //  return new TPZCompElPoint(mesh,gel,index);
}
static TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoLinear,TPZShapeLinear>(mesh,gel,index);
  //  return new TPZCompEl1d(mesh,gel,index);
}
static TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//  return new TPZIntelGen<TPZGeoQuad,TPZShapeQuad>(mesh,gel,index);
    return new TPZCompElQ2d(mesh,gel,index);
}
static TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>(mesh,gel,index);
  //  return new TPZCompElT2d(mesh,gel,index);
}
static TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoCube,TPZShapeCube>(mesh,gel,index);
  //  return new TPZCompElC3d(mesh,gel,index);
}
static TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoPrism,TPZShapePrism>(mesh,gel,index);
  //  return new TPZCompElPr3d(mesh,gel,index);
}
static TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoPyramid,TPZShapePiram>(mesh,gel,index);
  //  return new TPZCompElPi3d(mesh,gel,index);
}
static TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>(mesh,gel,index);
  //  return new TPZCompElT3d(mesh,gel,index);
}


TPZCompEl *(*TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePointEl;
TPZCompEl *(*TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;
TPZCompEl *(*TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateQuadEl;
TPZCompEl *(*TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;
TPZCompEl *(*TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateCubeEl;
TPZCompEl *(*TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePrismEl;
TPZCompEl *(*TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTetraEl;
TPZCompEl *(*TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePyramEl;
  */
