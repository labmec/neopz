/***************************************************************************
                          tpzgeoelrefpattern.cc  -  description
                             -------------------
    begin                : Tue Dec 23 2003
    copyright            : (C) 2003 by LabMeC - DES - FEC - UNICAMP (Edimar Cesar Rylo) & EMBRAER
    email                : cesar@labmec.fec.unicamp.br
 ***************************************************************************/

#include "tpzgeoelrefpattern.h"
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
#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern():TPZGeoElRefLess<TShape,TGeo>(){
  fSubEl = 0;
  fRefPattern = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::~TPZGeoElRefPattern(){
  if (fRefPattern) delete fRefPattern;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl = 0;
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh,index) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(0);
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(id,nodeindexes,matind,mesh) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(0);
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::Initialize : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index,TPZRefPattern *refpat) {
  TPZGeoElRefLess<TShape,TGeo>::Initialize(nodeindices,matind,mesh,index);
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(0);
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::Initialize : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::SetSubElement(int id, TPZGeoEl *el){
  int nsubel = fRefPattern->NSubElements();
  if (id<0 || id >nsubel){
    PZError << "TPZGeoElRefPattern::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id] = el;
  return;
}

template<class TShape, class TGeo>
REAL TPZGeoElRefPattern<TShape,TGeo>::RefElVolume(){
  return TShape::RefElVolume();
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::MidSideNodeIndex(int side,int &index){
//  TRef::MidSideNodeIndex(this,side,index);
}

template<class TShape, class TGeo>
int TPZGeoElRefPattern<TShape,TGeo>::NSubElements(){
  return fRefPattern->NSubElements();
}

template<class TShape, class TGeo>
int TPZGeoElRefPattern<TShape,TGeo>::NSideSubElements2(int side){
  return fRefPattern->NSideSubElements(side);
}

template<class TShape, class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TShape,TGeo>::SubElement(int is){
  int nsubel = fRefPattern->NSubElements();
  if(is<0 || is>nsubel){
    cout << "TPZGeoElRefPattern::SubElement index error is= " << is << endl;;
  }
  return fSubEl[is];
}

template<class TShape, class TGeo>
TPZGeoElSide TPZGeoElRefPattern<TShape,TGeo>::SideSubElement(int side,int position){
//  TPZStack<TPZGeoElSide> subs;
//  
//  TRef::GetSubElements(this,side,subs);
//  return subs[position];
  int sub, sideout;
  fRefPattern->SideSubElement(side,position,sub,sideout);
  if (!fSubEl[sub]) {
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::SideSubElement : Error subelement not found for side "
            << side << " position " << position << endl;
    return TPZGeoElSide();
  }
  return TPZGeoElSide (fSubEl[sub],sideout);
}

template<class TShape, class TGeo>
TPZTransform TPZGeoElRefPattern<TShape,TGeo>::GetTransform(int side,int son){
//  return TRef::GetTransform(side,son);
  return fRefPattern->Transform(side,son);  
}

template<class TShape, class TGeo>
void
TPZGeoElRefPattern<TShape,TGeo>::Divide(TPZVec<TPZGeoEl *> &pv){

//  TRef::Divide(this,pv);
}

template<class TShape, class TGeo>
void
TPZGeoElRefPattern<TShape,TGeo>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
  //TRef::GetSubElements(this,side,subel);
  int i,nsidesubel = fRefPattern->NSideSubElements(side);
  for (i=0;i<nsidesubel;i++){
    subel.Push(SideSubElement(side,i));  
  }
}

template class TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>;
template class TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>;
template class TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>;
template class TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>;
template class TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>;
template class TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>;
template class TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>;
template class TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>;
