// -*- c++ -*-

// $Id: pzelctemp.cpp,v 1.9 2003-11-06 14:36:10 cedric Exp $

#include "pzelctemp.h"
#include "pzquad.h"
#include "pzgeoel.h"

//template<class TGEO, class TSHAPE>
//class TPZIntelGen : public TPZInterpolatedElement {


template<class TGEO, class TSHAPE>
TPZIntelGen<TGEO,TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
  TPZInterpolatedElement(mesh,gel,index) {
  int i;
  for(i=0;i<TSHAPE::NSides-TSHAPE::NNodes;i++) {
    //    fSideOrder[i] = gOrder;
    fPreferredSideOrder = gOrder;
  }
  for(i=0; i<TSHAPE::NSides; i++) fConnectIndexes[i]=-1;
  //  RemoveSideRestraintsII(EInsert);
  gel->SetReference(this);
  for(i=0;i<TSHAPE::NNodes;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<TSHAPE::NSides;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
    IdentifySideOrder(i);
  }
  
  TPZManVector<int,3> order(3,2*TPZCompEl::gOrder);
  fIntRule.SetOrder(order);
  
}

template<class TGEO, class TSHAPE>
TPZIntelGen<TGEO,TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TGEO,TSHAPE> &copy) :
  TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
  fPreferredSideOrder = copy.fPreferredSideOrder;
  int i;
  for(i=0;i<TSHAPE::NSides;i++) {
    fConnectIndexes[i] = copy.fConnectIndexes[i];
  }
}

template<class TGEO, class TSHAPE>
TPZIntelGen<TGEO,TSHAPE>::~TPZIntelGen(){
  if(fReference) {
    if(fReference->Reference()) {
      RemoveSideRestraintsII(EDelete);
    }
    fReference->ResetReference();
  }
}

template<class TGEO, class TSHAPE>
MElementType TPZIntelGen<TGEO,TSHAPE>::Type() {
  return TGEO::Type();
}


//template<class TGEO, class TSHAPE>
//int TPZIntelGen<TGEO,TSHAPE>::NConnects() {
//  return TSHAPE::NSides;
//}

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::SetConnectIndex(int i, int connectindex){
#ifndef NDEBUG
  if(i<0 || i>= TSHAPE::NSides) {
    cout << " TPZIntelGen<TGEO,TSHAPE>::SetConnectIndex index " << i <<
      " out of range\n";
    return;
  }
#endif
  fConnectIndexes[i] = connectindex;
}

template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::NConnectShapeF(int connect){
  if(connect < TSHAPE::NNodes) return 1;
  int order = SideOrder(connect);
  if(order < 0) return 0;
  return TSHAPE::NConnectShapeF(connect, order);
}

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::SetIntegrationRule(int ord) {
  TPZManVector<int,3> order(TSHAPE::Dimension,ord);
  fIntRule.SetOrder(order);
}

//template<class TGEO, class TSHAPE>
//int TPZIntelGen<TGEO,TSHAPE>::Dimension() {
//    return TSHAPE::Dimension;
//}


//template<class TGEO, class TSHAPE>
//int TPZIntelGen<TGEO,TSHAPE>::NCornerConnects() {
//  return TSHAPE::NNodes;
//}

template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::NSideConnects(int side){
  return TSHAPE::NSideConnects(side);
}

template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::SideConnectLocId(int node, int side) {
  return TSHAPE::SideConnectLocId(side,node);
}

/**Sets the interpolation order for the interior of the element*/
template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::SetInterpolationOrder(int order) {
  fPreferredSideOrder = order;
}

/**Identifies the interpolation order on the interior of the element*/
template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(TSHAPE::NSides-TSHAPE::NNodes);
  int i;
  for(i=0; i<TSHAPE::NSides-TSHAPE::NNodes; i++) {
    ord[i] = SideOrder(i);
  }
}

/**return the preferred order of the polynomial along side iside*/
template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::PreferredSideOrder(int side) {
  if(side < TSHAPE::NNodes) return 0;
  if(side<TSHAPE::NSides) {
	  int order =fPreferredSideOrder;
	  return AdjustPreferredSideOrder(side,order);
  }
  PZError << "TPZIntelgen::PreferredSideOrder called for side = " << side << "\n";
  return 0;
  
}

template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::ConnectIndex(int con) {

#ifndef NDEBUG
  if(con<0 || con>= TSHAPE::NSides) {
    cout << "TPZIntelgen::ConnectIndex wrong parameter con " << con << 
      " NSides " << TSHAPE::NSides << endl;
    return -1;
  }

#endif
  return fConnectIndexes[con];
}



/**Sets the preferred interpolation order along a side
   This method only updates the datastructure of the element
   In order to change the interpolation order of an element, use the method PRefine
*/
template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::SetPreferredSideOrder(int order) {
  fPreferredSideOrder = order;
}

/**sets the interpolation order of side to order*/
template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::SetSideOrder(int side, int order) {
  if(side<0 || side >= TSHAPE::NSides || order <1) {
    PZError << "TPZIntelGen::SetSideOrder. Bad paramenter side.\n";
    return;
  }
  if(side>= TSHAPE::NNodes) {
    if(fConnectIndexes[side] == -1) return;
    TPZConnect &c = Connect(side);
    c.SetOrder(order);
    int seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial *mat = Material();
    if(mat) nvar = mat->NStateVariables();
    Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
    if(side == TSHAPE::NSides-1) {
      SetIntegrationRule(order);
    }
  }
}

/**returns the actual interpolation order of the polynomial along the side*/
template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::SideOrder(int side) {
  if(side < TSHAPE::NNodes || side >= TSHAPE::NSides) return 0;
  if(fConnectIndexes[side] == -1) return -1;
  TPZConnect &c = Connect(side);
  return c.Order();  
}

/**transform a point in the parameter space of the side into a point in the space
   of the master element*/
//template<class TGEO, class TSHAPE>
//void TPZIntelGen<TGEO,TSHAPE>::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {
  
//}

/**transform a point in the parameter space of the master element into a point in the
   space of the side*/
//virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);


/**compute the values of the shape function of the side*/

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {

}

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZManVector<int,TSHAPE::NNodes> id(TSHAPE::NNodes,0);
  TPZManVector<int, TSHAPE::NSides-TSHAPE::NNodes> ord(TSHAPE::NSides-TSHAPE::NNodes,0);
  int i;
  for(i=0; i<TSHAPE::NNodes; i++) {
    id[i] = fReference->NodePtr(i)->Id();
  }
  for(i=0; i<TSHAPE::NSides-TSHAPE::NNodes; i++) {
    ord[i] = SideOrder(i+TSHAPE::NNodes);
  }
  TSHAPE::Shape(pt,id,ord,phi,dphi);
}

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TSHAPE::Dimension && Material()->Id() > 0) {
    new typename TGEO::GraphElType(this,&grafgrid);
  }
}

// template<class TGEO, class TSHAPE>
// void TPZIntelGen<TGEO,TSHAPE>::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {

// }

  /** Jorge 09/06/2001
   * Returns the transformation which transform a point from the side to the interior of the element
   */
template<class TGEO, class TSHAPE>
TPZTransform TPZIntelGen<TGEO,TSHAPE>::TransformSideToElement(int side){
  return TSHAPE::TransformSideToElement(side);
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
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"

void TPZIntelGen<TPZGeoPoint,TPZShapePoint>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  phi(0,0) = 1.;
}

void TPZIntelGen<TPZGeoPoint,TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 0) cout << "A point element has no graphical representation\n";
}

void TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 3) cout << "A tetrahedra element has no graphical representation\n";
}

void TPZIntelGen<TPZGeoPrism,TPZShapePrism>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 3) cout << "A prism element has no graphical representation\n";
}

void TPZIntelGen<TPZGeoPyramid,TPZShapePiram>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 3) cout << "A pyramid element has no graphical representation\n";
}

template class TPZIntelGen<TPZGeoPoint,TPZShapePoint>;
template class TPZIntelGen<TPZGeoLinear, TPZShapeLinear>;
template class TPZIntelGen<TPZGeoQuad,TPZShapeQuad>;
template class TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>;
template class TPZIntelGen<TPZGeoCube,TPZShapeCube>;
template class TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>;
template class TPZIntelGen<TPZGeoPrism,TPZShapePrism>;
template class TPZIntelGen<TPZGeoPyramid,TPZShapePiram>;
