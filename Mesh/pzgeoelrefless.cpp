/***************************************************************************
                          pzgeoelrefless.cc  -  description
                             -------------------
    begin                : Fri Dec 12 2003
    copyright            : (C) 2003 by phil
    email                : phil@localhost
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "pzgeoelrefless.h"


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
//#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"
//#include "pzstack.h"

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess():TPZGeoEl(){
  int i;
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = -1;
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
//  fSubElement = -1;
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZGeoElRefLess  &gel):TPZGeoEl(gel){
  int i;
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = gel.fNodeIndexes[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i].SetConnectivity(gel.fNeighbours[i]);
//  fSubElement = -1;
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::~TPZGeoElRefLess(){
//  TPZGeoEl::~TPZGeoEl();
}
/** divides the element and puts the resulting elements in the vector */
//template<class TShape, class TGeo>
//void TPZGeoElRefLess<TShape,TGeo>::Divide(TPZVec < TPZGeoEl * > & pv){
//  pv.Resize (1);
//  TPZGeoEl *subel;
//  subel = new TPZGeoElRefLess(*this);
//  fSubElement = subel->Index();
//  subel->SetFather(this);
//  subel->SetFather(fIndex);
//}

/** return 1 if the element has subelements along side */
//template<class TShape, class TGeo>
//int TPZGeoElRefLess<TShape,TGeo>::HasSubElement(){
//  return 0;
//}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindices.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) :
  TPZGeoEl(matind,mesh,index) {

  int i,nnod = nodeindices.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(id,matind,mesh) {
  int i,nnod = nodeindexes.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template< class TShape, class TGeo>
void TPZGeoElRefLess<TShape,TGeo>::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index ) {
  int i;
  for(i = 0; i < TGeo::NNodes; i++ ){
     fNodeIndexes[ i ] = nodeindices[i];
  }
// for( int i = 0; i < TRef::NSubEl; i++ ){
//     fSubEl[ i ] = 0;
 // }
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NodeIndex(int node) {
  if(node<0 || node>7) return -1;
  return fNodeIndexes[node];
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideNodeIndex(int side,int node) {
  if(side<0 || side>(TShape::NSides - 1) || node<0) {
    PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return fNodeIndexes[TShape::SideNodeLocId(side,node)];
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideNodeLocIndex(int side,int node) {

  if(side<0 || side>(TShape::NSides - 1) || node<0) {
    PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return TShape::SideNodeLocId(side,node);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::SetSubElement(int id, TPZGeoEl *el){
  if (id!=0)
    PZError << "TPZGeoElRefLess<TShape,TGeo>::SetSubElement - Fodeu!\n";
  else PZError << "TPZGeoElRefLess<TShape,TGeo>::SetSubElement - Por enquanto eu no faco nada!\n";
  return;
}

template<class TShape, class TGeo>
TPZIntPoints *TPZGeoElRefLess<TShape,TGeo>::CreateSideIntegrationRule(int side, int order){
  return TGeo::CreateSideIntegrationRule(side,order);
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NNodes() {
  return TGeo::NNodes;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NCornerNodes(){
  return TGeo::NNodes;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSides(){
  return TGeo::NSides;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideNodeLocId(int side, int node){
  return TShape::SideNodeLocId(side,node);
}

template<class TShape, class TGeo>
REAL
TPZGeoElRefLess<TShape,TGeo>::RefElVolume(){
  return TShape::RefElVolume();
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSideNodes(int side){
  return TShape::NSideNodes(side);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::MidSideNodeIndex(int side,int &index){
  //TRef::MidSideNodeIndex(this,side,index);
  index = -1;
  if(side<0 || side>NSides()-1) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::MidSideNodeIndex. Bad parameter side = " << side << endl;
    return;
  }
  if(side<NNodes()) {//o nó medio do lado 0 é o 0 etc.
    index = NodeIndex(side);
    return;
  }
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideIsUndefined(int side){
  if (side < 0 || side > NSides()){
    PZError << "TPZGeoElRefLess<TShape,TGeo>::SideIsUndefined - bad side: " << side << endl;
  }
  return (fNeighbours[side].Side() == -1);
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSubElements(){
  //return TRef::NSubEl;
  return 0;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSideSubElements(int side){
  //return TRef::NSideSubElements(side);
  return 0;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSideSubElements2(int side){
  // return TRef::NSideSubElements(side);
  return 0;   
}

template<class TShape, class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TShape,TGeo>::CreateBCGeoEl(int side, int bc){
  return TGeo::CreateBCGeoEl(this,side,bc);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::SetNodeIndex(int i,int nodeindex){
  if(i<0 || i>(TGeo::NNodes - 1)){
    cout << "TPZGeoElRefLess::SetNodeIndex index error i = " << i << endl;
    return;
  }
  fNodeIndexes[i] = nodeindex;
}

template<class TShape, class TGeo>
TPZTransform
TPZGeoElRefLess<TShape,TGeo>::SideToSideTransform(int sidefrom,int sideto){
  return TShape::SideToSideTransform(sidefrom,sideto);
}

template<class TShape, class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TShape,TGeo>::SubElement(int is){
  if(is<0 || is>1){//(TRef::NSubEl - 1)){
    cout << "TPZGeoElRefLess::SubElement index error is= " << is << endl;;
  }
//  return fSubEl[is];
  return 0;
}
/* 
template<class TShape, class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TShape,TGeo>::SideSubElement(int side,int position){
 TPZStack<TPZGeoElSide> subs;
  TRef::GetSubElements(this,side,subs);
  return subs[position];
  if (side < 0 || side > NSides()){
    PZError << "TPZGeoElRefLess<TShape,TGeo>::SideSubElement - bad side: " << side << endl;
    return TPZGeoElSide();
  }
  return (fNeighbours[side].Side() == -1);
}  */

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideDimension(int side){
  return TShape::SideDimension(side);
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::Dimension(){
  return TShape::Dimension;
}

template<class TShape, class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TShape,TGeo>::HigherDimensionSides(int side,int targetdimension){
  cout << "TPZGeoElRefLess::HigherDimensionSides nao deve ser usado\n";
  return TPZGeoElSide();
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
  TPZStack<int> highsides;
  TShape::HigherDimensionSides(side,highsides);
  int i,size = highsides.NElements();
  for(i=0;i<size;i++) {
    if(SideDimension(highsides[i]) == targetdimension) {
      elsides.Push(TPZGeoElSide(this,highsides[i]));
    }
  }
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::LowerDimensionSides(int side,TPZStack<int> &smallsides){
  int nsidecon = TShape::NSideConnects(side);
  int is;
  for(is=0; is<nsidecon-1; is++)
    smallsides.Push(TShape::SideConnectLocId(side,is));
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::BuildTransform(int side, TPZGeoEl *father,TPZTransform &t){
  BuildTransform2(side,father,t);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix &jac,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
  TPZFNMatrix<3*TGeo::NNodes> nodes(3,TGeo::NNodes);
  TPZGeoNode *np;
  TPZAdmChunkVector<TPZGeoNode> &nodevec = Mesh()->NodeVec();
  int i,j;
  for(i=0;i<TGeo::NNodes;i++) {
    np = &nodevec[fNodeIndexes[i]];
    for(j=0;j<3;j++) {
      nodes(j,i) = np->Coord(j);
    }
  }
  TGeo::Jacobian(nodes,coordinate,jac,axes,detjac,jacinv);
//   if(TGeo::NNodes == 2) {
//     detjac = 1.;
//     jacinv(0,0) = 1.;
//     jac(0,0) = 1.;
//   }
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result){
  TPZFNMatrix<3*TGeo::NNodes> nodes(3,TGeo::NNodes);
  TPZGeoNode *np;
  TPZAdmChunkVector<TPZGeoNode> &nodevec = Mesh()->NodeVec();
  int i,j;
  for(i=0;i<TGeo::NNodes;i++) {
    np = &nodevec[fNodeIndexes[i]];
    for(j=0;j<3;j++) {
      nodes(j,i) = np->Coord(j);
    }
  }
  TGeo::X(nodes,coordinate,result);
}

template<class TShape, class TGeo>
TPZTransform
TPZGeoElRefLess<TShape,TGeo>::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t)
{
  //Augusto:09/01/01
  TPZGeoEl *myfather = Father();
  if(side<0 || side>(TShape::NSides-1) || !myfather){
    PZError << "TPZGeoElRefLess::BuildTransform2 side out of range or father null\n";
    return TPZTransform(0,0);
  }
  TPZGeoElSide fathloc = Father2(side);
  int son = WhichSubel();
  TPZTransform trans=myfather->GetTransform(side,son);
  trans = trans.Multiply(t);
  if(fathloc.Element() == father) return trans;
  trans = myfather->BuildTransform2(fathloc.Side(),father,trans);
  return trans;
}


template<class TShape, class TGeo>
TPZTransform
TPZGeoElRefLess<TShape,TGeo>::GetTransform(int side,int son){
//  return TRef::GetTransform(side,son);
//  if(side<0 || side>NSides()-1){
    PZError << "TPZGeoElRefLess<TShape,TGeo>::GetTransform::Never should be called\n";
    return TPZTransform(0,0);
//  }
/*  int smalldim = TShape::SideDimension(side);
  int fatherside = FatherSide(side,son);
  int largedim = TShape::SideDimension(fatherside);
  TPZTransform trans(largedim,smalldim);
  int i,j;
  for(i=0; i<largedim; i++) {
    for(j=0; j<smalldim; j++) {
      trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
    }
    trans.Sum() (i,0) = buildt[whichsubel][side][2][i];
  }
  return trans;
  */
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::CenterPoint(int side, TPZVec<REAL> &cent){
  TShape::CenterPoint(side,cent);
}

template<class TShape, class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TShape,TGeo>::Father2(int side)
{
  //cout << " Father2 teste Cedric: 08/05/2003\n";
  TPZGeoEl *father = Father();
  if(!father) return TPZGeoElSide();
  int son = WhichSubel();
  if(son<0) return TPZGeoElSide();
  int fathsid = father->FatherSide(side,son);
  return TPZGeoElSide(father,fathsid);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
//  TRef::GetSubElements(this,side,subel);
  return;
}

template class TPZGeoElRefLess<TPZShapeCube,TPZGeoCube>;
template class TPZGeoElRefLess<TPZShapeLinear,TPZGeoLinear>;
template class TPZGeoElRefLess<TPZShapeQuad,TPZGeoQuad>;
template class TPZGeoElRefLess<TPZShapeTriang,TPZGeoTriangle>;
template class TPZGeoElRefLess<TPZShapePrism,TPZGeoPrism>;
template class TPZGeoElRefLess<TPZShapeTetra,TPZGeoTetrahedra>;
template class TPZGeoElRefLess<TPZShapePiram,TPZGeoPyramid>;
template class TPZGeoElRefLess<TPZShapePoint,TPZGeoPoint>;

int main_refless(){
//TPZGeoElRefLess <TPZShapeCube,TPZGeoCube> el1;
//TPZGeoElRefLess <TPZShapeLinear,TPZGeoLinear> el2;
//TPZGeoElRefLess <TPZShapeQuad,TPZGeoQuad> el3;
//TPZGeoElRefLess <TPZShapeTriang,TPZGeoTriangle> el4;
//TPZGeoElRefLess <TPZShapePrism,TPZGeoPrism> el5;
//TPZGeoElRefLess <TPZShapeTetra,TPZGeoTetrahedra> el6;
//TPZGeoElRefLess <TPZShapePiram,TPZGeoPyramid> el7;
//TPZGeoElRefLess <TPZShapePoint,TPZGeoPoint> el8;

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
  return new TPZIntelGen<TPZGeoQuad,TPZShapeQuad>(mesh,gel,index);
//    return new TPZCompElQ2d(mesh,gel,index);
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


TPZCompEl *(*TPZGeoElRefLess<TPZShapePoint,TPZGeoPoint>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePointEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapeLinear,TPZGeoLinear>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapeQuad,TPZGeoQuad>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateQuadEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapeTriang,TPZGeoTriangle>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapeCube,TPZGeoCube>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateCubeEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapePrism,TPZGeoPrism>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePrismEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapeTetra,TPZGeoTetrahedra>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTetraEl;
TPZCompEl *(*TPZGeoElRefLess<TPZShapePiram,TPZGeoPyramid>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePyramEl;


