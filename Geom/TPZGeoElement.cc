#include "TPZGeoElement.h"
//#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"

template<class TShape, class TGeo, class TRef>
TPZCompEl *(*TPZGeoElement<TShape,TGeo,TRef>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);


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
	return -1;//TGeo::NSubElements();
}

template<class TShape, class TGeo, class TRef>
int
TPZGeoElement<TShape,TGeo,TRef>::NSideSubElements(int side){
	return -1;//TRef::NSideSubElements(side);
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
TPZGeoElement<TShape,TGeo,TRef>::AllHigherDimensionSides(int side,
			     int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> highsides;
	TShape::HigherDimensionSides(side,highsides);
	int i,size = highsides.NElements();
	for(i=0;i<size;i++) elsides.Push(TPZGeoElSide(this,highsides[i]));
}

template<class TShape, class TGeo, class TRef>
void
TPZGeoElement<TShape,TGeo,TRef>::LowerDimensionSides(int side,TPZStack<int> &smallsides){
	int nsidecon = TShape::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TShape::SideConnectLocId(side,is));
	return;

	//TPZStack<int> smallsidesint;
	//TShape::LowerDimensionSides(side,smallsidesint);
	//int i,size = smallsidesint.NElements();
	//for(i=0;i<size;i++) smallsides.Push(TPZGeoElSide(this,smallsidesint[i]));
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

	return TShape::CenterPoint(side,cent);
}

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


template class TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>;
template class TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>;
template class TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>;
template class TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>;
template class TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>;
template class TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>;
template class TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>;
template class TPZGeoElement<TPZShapeLinear,TPZGeoPoint,TPZRefPoint>;

int newteste(){
TPZGeoElement <TPZShapeCube,TPZGeoCube,TPZRefCube> el1;
TPZGeoElement <TPZShapeLinear,TPZGeoLinear,TPZRefLinear> el2;
TPZGeoElement <TPZShapeQuad,TPZGeoQuad,TPZRefQuad> el3;
TPZGeoElement <TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle> el4;
TPZGeoElement <TPZShapePrism,TPZGeoPrism,TPZRefPrism> el5;
TPZGeoElement <TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra> el6;
TPZGeoElement <TPZShapePiram,TPZGeoPyramid,TPZRefPyramid> el7;
TPZGeoElement <TPZShapeLinear,TPZGeoPoint,TPZRefPoint> el8;

return 0;
}
