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
  TPZGeoEl(matind,mesh) {
  
  int i,nnod = nodeindices.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElement<TShape,TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  for(i=0;i<TGeo::NNodes;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo, class TRef>
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) : 
  TPZGeoEl(matind,mesh,index) {
  
  int i,nnod = nodeindices.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElement<TShape,TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo, class TRef>
TPZGeoElement<TShape,TGeo,TRef>::TPZGeoElement(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(id,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElement<TShape,TGeo,TRef>::Constuctor, number of nodes : " << nnod << endl;
    return;
  }
  
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<TRef::NSubEl;i++) fSubEl[i] = 0;
}

/*
template< class TShape, class TGeo, class TRef >
TPZGeoElement< TShape, TGeo, TRef >::TPZGeoElement(
   int* nodeindices, int matind, TPZGeoMesh& mesh ) :
   TPZGeoEl( matind, mesh )
{
   // BEWARE! We do not test the number of elements here, cause we are
   // passing a pointer to some memory portion!
  
  for( int i = 0; i < TGeo::NNodes; i++ )
  {
     fNodeIndexes[ i ] = *nodeindices++ - 1;
  }

  for( int i = 0; i < TGeo::NNodes; i++ )
  {
     fSubEl[ i ] = 0;
  }
}
*/

/*
template< class TShape, class TGeo, class TRef >
TPZGeoElement< TShape, TGeo, TRef >::TPZGeoElement(
   int* nodeindices, int matind, TPZGeoMesh& mesh, int& index ) :
   TPZGeoEl( matind, mesh, index )
{
   // BEWARE! We do not test the number of elements here, cause we are
   // passing a pointer to some memory portion!
  
  for( int i = 0; i < TGeo::NNodes; i++ )
  {
     fNodeIndexes[ i ] = *nodeindices++ - 1;
  }

  for( int i = 0; i < TRef::NSubEl; i++ )
  {
     fSubEl[ i ] = 0;
  }
}
*/
template< class TShape, class TGeo, class TRef >
void TPZGeoElement< TShape, TGeo, TRef >::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index ) {
  
  for( int i = 0; i < TGeo::NNodes; i++ )
  {
     fNodeIndexes[ i ] = nodeindices[i];
  }

  for( int i = 0; i < TRef::NSubEl; i++ )
  {
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
int 
TPZGeoElement<TShape,TGeo,TRef>::NodeIndex(int node) {

  if(node<0 || node>7) return -1;
  return fNodeIndexes[node];
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::SideNodeIndex(int side,int node) {

  if(side<0 || side>(TShape::NSides - 1) || node<0) {
    PZError << "TPZGeoElement::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return fNodeIndexes[TShape::SideNodeLocId(side,node)];
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::SideNodeLocIndex(int side,int node) {

  if(side<0 || side>(TShape::NSides - 1) || node<0) {
    PZError << "TPZGeoElement::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return TShape::SideNodeLocId(side,node);
}

template<class TShape, class TGeo, class TRef>
void 
TPZGeoElement<TShape,TGeo,TRef>::SetSubElement(int id, TPZGeoEl *el){

  if (id<0 || id >(TRef::NSubEl - 1)){
    PZError << "TPZGeoElement::Trying do define subelement :" 
	    << id << endl;
    return;
  }
  fSubEl[id] = el;
  return;
}

template<class TShape, class TGeo, class TRef>
TPZIntPoints *TPZGeoElement<TShape,TGeo,TRef>::CreateSideIntegrationRule(int side, int order){
	return TGeo::CreateSideIntegrationRule(side,order);
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NNodes() {
	return TGeo::NNodes;
}

template<class TShape, class TGeo, class TRef>
int 
TPZGeoElement<TShape,TGeo,TRef>::NCornerNodes(){
	return TGeo::NNodes;
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSides(){
	return TGeo::NSides;
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::SideNodeLocId(int side, int node){
	return TShape::SideNodeLocId(side,node);
}

template<class TShape, class TGeo, class TRef>
REAL
TPZGeoElement<TShape,TGeo,TRef>::RefElVolume(){
	return TShape::RefElVolume();
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSideNodes(int side){
	return TShape::NSideNodes(side);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::MidSideNodeIndex(int side,int &index){
	TRef::MidSideNodeIndex(this,side,index);
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::SideIsUndefined(int side){
	return (fNeighbours[side].Side() == -1);
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSubElements(){
	return TRef::NSubEl;
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSideSubElements(int side){
	return TRef::NSideSubElements(side);
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSideSubElements2(int side){
	return TRef::NSideSubElements(side);
}

template<class TShape, class TGeo, class TRef>
TPZGeoEl *
TPZGeoElement<TShape,TGeo,TRef>::CreateBCGeoEl(int side, int bc){
  return TGeo::CreateBCGeoEl(this,side,bc);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::SetNodeIndex(int i,int nodeindex){

  if(i<0 || i>(TGeo::NNodes - 1)){
    cout << "TPZGeoElement::SetNodeIndex index error i = " << i << endl;
    return;
  }
  fNodeIndexes[i] = nodeindex;
}

template<class TShape, class TGeo, class TRef>
TPZTransform
TPZGeoElement<TShape,TGeo,TRef>::SideToSideTransform(int sidefrom,int sideto){
  return TShape::SideToSideTransform(sidefrom,sideto);
}

template<class TShape, class TGeo, class TRef>
TPZGeoEl *
TPZGeoElement<TShape,TGeo,TRef>::SubElement(int is){
  if(is<0 || is>(TRef::NSubEl - 1)){
    cout << "TPZGeoElement::SubElement index error is= " << is << endl;;
  }
  return fSubEl[is];
}

template<class TShape, class TGeo, class TRef>
TPZGeoElSide
TPZGeoElement<TShape,TGeo,TRef>::SideSubElement(int side,int position){
  TPZStack<TPZGeoElSide> subs;
  TRef::GetSubElements(this,side,subs);
  return subs[position];
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::SideDimension(int side){
  return TShape::SideDimension(side);
}

template<class TShape, class TGeo, class TRef>  
int
TPZGeoElement<TShape,TGeo,TRef>::Dimension(){
  return TShape::Dimension;
}

template<class TShape, class TGeo, class TRef> 
TPZGeoElSide
TPZGeoElement<TShape,TGeo,TRef>::HigherDimensionSides(int side,int targetdimension){
	cout << "TPZGeoElement::HigherDimensionSides nao deve ser usado\n";
	return TPZGeoElSide();
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
  TPZStack<int> highsides;
  TShape::HigherDimensionSides(side,highsides);
  int i,size = highsides.NElements();
  for(i=0;i<size;i++)
    elsides.Push(TPZGeoElSide(this,highsides[i]));
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::LowerDimensionSides(int side,TPZStack<int> &smallsides){
	int nsidecon = TShape::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++)
	  smallsides.Push(TShape::SideConnectLocId(side,is));
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::BuildTransform(int side, TPZGeoEl *father,TPZTransform &t){
	BuildTransform2(side,father,t);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix &jac,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
  TPZFMatrix nodes(3,TGeo::NNodes);
  TPZGeoNode *np;
  int i,j;
  for(i=0;i<TGeo::NNodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++) {
      nodes(j,i) = np->Coord(j);
    }
  }
  TGeo::Jacobian(nodes,coordinate,jac,axes,detjac,jacinv);	
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result){
  TPZFMatrix nodes(3,TGeo::NNodes);
  TPZGeoNode *np;
  int i,j;
  for(i=0;i<TGeo::NNodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++) {
      nodes(j,i) = np->Coord(j);
    }
  }
  TGeo::X(nodes,coordinate,result);
}

template<class TShape, class TGeo, class TRef>
TPZTransform
TPZGeoElement<TShape,TGeo,TRef>::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t){//Augusto:09/01/01
  
  if(side<0 || side>(TShape::NSides-1) || !fFather){
    PZError << "TPZGeoElement::BuildTransform2 side out of range or father null\n";
    return TPZTransform(0,0);
  }
  TPZGeoElSide fathloc = Father2(side);
  int son = WhichSubel();
  TPZTransform trans=fFather->GetTransform(side,son);
  trans = trans.Multiply(t);
  if(fathloc.Element() == father) return trans;
  trans = fFather->BuildTransform2(fathloc.Side(),father,trans);
  return trans;
}


template<class TShape, class TGeo, class TRef>
TPZTransform
TPZGeoElement<TShape,TGeo,TRef>::GetTransform(int side,int son){

	return TRef::GetTransform(side,son);
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::CenterPoint(int side, TPZVec<REAL> &cent){

  TShape::CenterPoint(side,cent);
}

template<class TShape, class TGeo, class TRef>
TPZGeoElSide 
TPZGeoElement<TShape,TGeo,TRef>::Father2(int side){//cout << " Father2 teste Cedric: 08/05/2003\n";
  if(!fFather) return TPZGeoElSide();
  int son = WhichSubel(); 
  if(son<0) return TPZGeoElSide();
  int fathsid = fFather->FatherSide(side,son);
  return TPZGeoElSide(fFather,fathsid);     
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
  return new TPZIntelGen<TPZGeoQuad,TPZShapeQuad>(mesh,gel,index);
  //  return new TPZCompElQ2d(mesh,gel,index);
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
