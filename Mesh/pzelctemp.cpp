// -*- c++ -*-

// $Id: pzelctemp.cpp,v 1.32 2007-04-12 14:16:59 cesar Exp $

#include "pzelctemp.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzintelgen"));
#endif

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

  int sideorder = SideOrder(TSHAPE::NSides-1);
  sideorder = 2*sideorder;
  if (sideorder > fIntRule.GetMaxOrder()) sideorder = fIntRule.GetMaxOrder();
  //  TPZManVector<int,3> order(3,2*sideorder+2);
  TPZManVector<int,3> order(3,sideorder);
  //TPZManVector<int,3> order(3,20);
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
cerviTPZIntelGen<TGEO,TSHAPE>::TPZIntelGen(TPZCompMesh &mesh,
                                      const TPZIntelGen<TGEO,TSHAPE> &copy,
                                      std::map<int,int> & gl2lcConMap,
                                      std::map<int,int> & gl2lcElMap) :
    TPZInterpolatedElement(mesh,copy,gl2lcElMap), fIntRule(copy.fIntRule)
{
//   if (gl2lcElMap.find(copy.fIndex)!= gl2lcElMap.end())
//   {
//     std::stringstream sout;
//     sout << "ERROR in : " << __PRETTY_FUNCTION__
//         << " trying to clone an already cloned element index: " << copy.fIndex;
//     LOGPZ_ERROR(logger, sout.str().c_str());
//     exit(-1);
//   }

  fPreferredSideOrder = copy.fPreferredSideOrder;
  int i;
  for(i=0;i<TSHAPE::NSides;i++)
  {
    int lcIdx = -1;
    int glIdx = copy.fConnectIndexes[i];
    if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
    else
    {
      std::stringstream sout;
      sout << "ERROR in : " << __PRETTY_FUNCTION__
          << " trying to clone the connect index: " << glIdx
          << " wich is not in mapped connect indexes!";
      LOGPZ_ERROR(logger, sout.str().c_str());
      fConnectIndexes[i] = -1;
      return;
    }
    fConnectIndexes[i] = lcIdx;
  }
//   gl2lcElMap[copy.fIndex] = this->Index();
}


template<class TGEO, class TSHAPE>
TPZIntelGen<TGEO,TSHAPE>::TPZIntelGen() :
  TPZInterpolatedElement(), fIntRule() {
  fPreferredSideOrder = -1;
  int i;
  for(i=0;i<TSHAPE::NSides;i++) {
    fConnectIndexes[i] = -1;
  }
}

template<class TGEO, class TSHAPE>
TPZIntelGen<TGEO,TSHAPE>::~TPZIntelGen(){
  if(Reference()) {
    if(Reference()->Reference()) {
      RemoveSideRestraintsII(EDelete);
    }
    Reference()->ResetReference();
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
#ifndef NODEBUG
  if(i<0 || i>= TSHAPE::NSides) {
    std::cout << " TPZIntelGen<TGEO,TSHAPE>::SetConnectIndex index " << i <<
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
    ord[i] = SideOrder(i+TSHAPE::NNodes);
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

#ifndef NODEBUG
  if(con<0 || con>= TSHAPE::NSides) {
    std::cout << "TPZIntelgen::ConnectIndex wrong parameter con " << con <<
      " NSides " << TSHAPE::NSides << std::endl;
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
  if(side<0 || side >= TSHAPE::NSides || (side >= TSHAPE::NNodes && order <1)) {
    PZError << "TPZIntelGen::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
    return;
  }
  if(side>= TSHAPE::NNodes) {
    if(fConnectIndexes[side] == -1) return;
    TPZConnect &c = Connect(side);
    c.SetOrder(order);
    int seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZAutoPointer<TPZMaterial> mat = Material();
    if(mat) nvar = mat->NStateVariables();
    Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
    if(side == TSHAPE::NSides-1) {
      SetIntegrationRule(2*order);
    }
  }
}

/**returns the actual interpolation order of the polynomial along the side*/
template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::SideOrder(int side) {
  if(side < TSHAPE::NNodes || side >= TSHAPE::NSides) return 0;
  if(fConnectIndexes[side] == -1)
  {
    std::cout << __PRETTY_FUNCTION__ << " side " << side << std::endl;
    Print(cout);
    return -1;
  }
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

  int nc = TSHAPE::NSideConnects(side);
  int nn = TSHAPE::NSideNodes(side);
  TPZManVector<int,27> id(nn),order(nc-nn);
  int n,c;
  TPZGeoEl *ref = Reference();
  for (n=0;n<nn;n++){
    int nodloc = TSHAPE::SideNodeLocId(side,n);
    id [n] = ref->NodePtr(nodloc)->Id();
  }
  for (c=nn;c<nc;c++){
    int conloc = TSHAPE::SideConnectLocId(side,c);
    order[c-nn] = SideOrder(conloc);
  }
  TSHAPE::SideShape(side, point, id, order, phi, dphi);
}

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZManVector<int,TSHAPE::NNodes> id(TSHAPE::NNodes,0);
  TPZManVector<int, TSHAPE::NSides-TSHAPE::NNodes+1> ord(TSHAPE::NSides-TSHAPE::NNodes,0);
  int i;
  TPZGeoEl *ref = Reference();
  for(i=0; i<TSHAPE::NNodes; i++) {
    id[i] = ref->NodePtr(i)->Id();
  }
  for(i=0; i<TSHAPE::NSides-TSHAPE::NNodes; i++) {
    ord[i] = SideOrder(i+TSHAPE::NNodes);
  }
  TSHAPE::Shape(pt,id,ord,phi,dphi);
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

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
/*
template<class TGEO, class TSHAPE>
int TPZIntelGen<TGEO,TSHAPE>::ClassId() const
{
  std::cout << "TPZIntelGen<TGEO,TSHAPE>::ClassId() type not specified at " << __FILE__ << ':' << __LINE__ << std::endl;
  return -1;
}
*/
  /**
  Save the element data to a stream
  */
template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::Write(TPZStream &buf, int withclassid)
{
  TPZInterpolatedElement::Write(buf,withclassid);
  TPZManVector<int,3> order(3,0);
  fIntRule.GetOrder(order);
  WriteObjects(buf,order);
  buf.Write(fConnectIndexes,TSHAPE::NSides);
  buf.Write(&fPreferredSideOrder,1);
}

  /**
  Read the element data from a stream
  */
template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::Read(TPZStream &buf, void *context)
{
  TPZInterpolatedElement::Read(buf,context);
  TPZManVector<int,3> order;
  ReadObjects(buf,order);
  fIntRule.SetOrder(order);
  buf.Read(fConnectIndexes,TSHAPE::NSides);
  buf.Read(&fPreferredSideOrder,1);
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
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"

/*
void TPZIntelGen<TPZGeoPoint,TPZShapePoint>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  phi(0,0) = 1.;
}
*/
using namespace pzgeom;
using namespace pzshape;

template<>
void TPZIntelGen<TPZGeoPoint,TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}

template<>
void TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 3)
  {
    new TPZGeoTetrahedra::GraphElType(this,&grafgrid);
  }
}

template<>
void TPZIntelGen<TPZGeoPrism,TPZShapePrism>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 3) std::cout << "A prism element has no graphical representation\n";
}

template<>
void TPZIntelGen<TPZGeoPyramid,TPZShapePiram>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 3) std::cout << "A pyramid element has no graphical representation\n";
}
#ifndef DOS

template<class TGEO, class TSHAPE>
void TPZIntelGen<TGEO,TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TSHAPE::Dimension && Material()->Id() > 0) {
    new typename TGEO::GraphElType(this,&grafgrid);
  }
}

#else
template<>
void TPZIntelGen<TPZGeoLinear,TPZShapeLinear>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TPZShapeLinear::Dimension && Material()->Id() > 0) {
    new typename TPZGeoLinear::GraphElType(this,&grafgrid);
  }
}
template<>
void TPZIntelGen<TPZGeoQuad,TPZShapeQuad>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TPZShapeQuad::Dimension && Material()->Id() > 0) {
    new typename TPZGeoQuad::GraphElType(this,&grafgrid);
  }
}
template<>
void TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TPZShapeTriang::Dimension && Material()->Id() > 0) {
    new typename TPZGeoTriangle::GraphElType(this,&grafgrid);
  }
}
template<>
void TPZIntelGen<TPZGeoCube,TPZShapeCube>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TPZShapeCube::Dimension && Material()->Id() > 0) {
    new typename TPZGeoCube::GraphElType(this,&grafgrid);
  }
}

#endif
template<>
int TPZIntelGen<TPZGeoPoint,TPZShapePoint>::ClassId() const
{
  return TPZINTELPOINTID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoPoint,TPZShapePoint>, TPZINTELPOINTID>;


template<>
int TPZIntelGen<TPZGeoLinear,TPZShapeLinear>::ClassId() const
{
  return TPZINTELLINEARID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoLinear,TPZShapeLinear>, TPZINTELLINEARID>;

template<>
int TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>::ClassId() const
{
  return TPZINTELTRIANGLEID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>, TPZINTELTRIANGLEID>;

template<>
int TPZIntelGen<TPZGeoQuad,TPZShapeQuad>::ClassId() const
{
  return TPZINTELQUADID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoQuad,TPZShapeQuad>, TPZINTELQUADID>;

template<>
int TPZIntelGen<TPZGeoCube,TPZShapeCube>::ClassId() const
{
  return TPZINTELCUBEID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoCube,TPZShapeCube>, TPZINTELCUBEID>;

template<>
int TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>::ClassId() const
{
  return TPZINTELTETRAID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>, TPZINTELTETRAID>;

template<>
int TPZIntelGen<TPZGeoPrism,TPZShapePrism>::ClassId() const
{
  return TPZINTELPRISMID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoPrism,TPZShapePrism>, TPZINTELPRISMID>;

template<>
int TPZIntelGen<TPZGeoPyramid,TPZShapePiram>::ClassId() const
{
  return TPZINTELPYRAMID;
}
template class
    TPZRestoreClass< TPZIntelGen<TPZGeoPyramid,TPZShapePiram>, TPZINTELPYRAMID>;


template class TPZIntelGen<TPZGeoPoint,TPZShapePoint>;
template class TPZIntelGen<TPZGeoLinear, TPZShapeLinear>;
template class TPZIntelGen<TPZGeoQuad,TPZShapeQuad>;
template class TPZIntelGen<TPZGeoTriangle,TPZShapeTriang>;
template class TPZIntelGen<TPZGeoCube,TPZShapeCube>;
template class TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra>;
template class TPZIntelGen<TPZGeoPrism,TPZShapePrism>;
template class TPZIntelGen<TPZGeoPyramid,TPZShapePiram>;
