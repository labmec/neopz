// -*- c++ -*-
#include "pzelcpoint.h"
#include "pzelct2d.h"
#include "pzelcq2d.h"
#include "pzerror.h"
#include "pzelgt2d.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzconnect.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzcmesh.h"
#include "pzelgpoint.h"
#include "pztrigraph.h"
#include "pzgraphel.h"



TPZCompElPoint::TPZCompElPoint(TPZCompMesh &mesh,TPZGeoElPoint *ref,int &index) :
		TPZInterpolatedElement(mesh,ref,index), fIntRule() {
  fConnectIndexes = -1;
  // Jorge 19/5/99
  /** Compatibilizing the continuity over the neighboard elements */
  TPZStack<TPZCompElSide> elvec;
  TPZCompElSide thisside(this,0);
  thisside.EqualLevelElementList(elvec,0,0);
  /** Asking if exist neighbour : If this neighboard is continuous all neighboards are continuous */
  if(elvec.NElements()) {
    ((TPZInterpolatedElement *)elvec[0].Element())->MakeConnectContinuous(elvec[0].Side());
  }

//  RemoveSideRestraintsII(EInsert);
  ref->SetReference(this);
  fConnectIndexes = CreateMidSideConnect(0);
  mesh.ConnectVec()[fConnectIndexes].IncrementElConnected();
}

TPZCompElPoint::TPZCompElPoint(TPZCompMesh &mesh,TPZGeoElPoint *ref,int &index,int /*noconnects*/) :
		TPZInterpolatedElement(mesh,ref,index), fIntRule() {

  fConnectIndexes = -1;
//  RemoveSideRestraintsII(EInsert);
}

void TPZCompElPoint::SetConnectIndex(int i,int connectindex) {
  if(i == 0 ) fConnectIndexes = connectindex;
  else {
    PZError << "TPZCompElPoint::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElPoint::ConnectIndex(int i) {
  if(i != 0) {
    PZError << "TPZCompElPoint::ConnectIndex. Bad parameter i.\n";
    return -1;
  }
  return fConnectIndexes;
}

int TPZCompElPoint::NConnectShapeF(int side) {
	  if(side==0)  return 1;
     PZError << "TPZCompElPoint::NConnectShapeF, bad parameter iconnect " << side << endl;
     return 0;
}

int TPZCompElPoint::NSideConnects(int side) {
  if(side!=0) {
    PZError << "TPZCompElPoint::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  return 1;
}

/**It do not verify the values of the c*/
int TPZCompElPoint::SideConnectLocId(int c, int side) {

    if(c==0 && side==0) return 0;
    PZError << "TPZCompElPoint::SideConnectLocId, connect = " << c << endl;
    return -1;
}

void TPZCompElPoint::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=1)
    PZError << "TPZCompElPoint::SetInterpolationOrder. ord has bad number of elements.\n";
}

void TPZCompElPoint::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(1);
  ord[0] = -1;
}

TPZIntPoints *TPZCompElPoint::CreateSideIntegrationRule(int ) {
   return new TPZInt1Point();
}

int TPZCompElPoint::PreferredSideOrder(int /*side*/) {
  return 0;
}

void TPZCompElPoint::SetPreferredSideOrder(int side, int /*order*/) {
  if (side)
    PZError << "TPZCompElPoint::SetPreferredSideOrder. Bad paramenter side: " << side << "\n";
  //  if(side != 0) return;
}


void TPZCompElPoint::SetSideOrder(int /*side*/, int order) {
  if(order) PZError << "TPZCompElPoint::SetSideOrder. Bad paramenter side.\n";
}

int TPZCompElPoint::SideOrder(int /*side*/) {
  return 0;
}

void TPZCompElPoint::SideParameterToElement(int side,TPZVec<REAL> &,TPZVec<REAL> &/*point*/) {
  if(side==0) return;
  PZError << "TPZCompElPoint::SideParameterToElement. Bad paramenter side.\n";
}

void TPZCompElPoint::ElementToSideParameter(int side,TPZVec<REAL> &,TPZVec<REAL> &) {

  if(side==0) return;
  PZError << "TPZCompElPoint::ElementToSideParameter. Bad side = " << side << ".\n";
}

void TPZCompElPoint::CornerShape(TPZVec<REAL> &/*pt*/,TPZFMatrix &phi,TPZFMatrix &dphi) {
  phi.Redim(1,1);
  phi(0,0) = 1.;
  dphi.Redim(0,1);
}

void TPZCompElPoint::SideShapeFunction(int side,TPZVec<REAL> &,TPZFMatrix &phi,TPZFMatrix &dphi) {  //    id ???
  if(side!=0) {
    PZError << "TPZCompElPoint::SideShapeFunction. Bad parameter side.\n";
    return;
  }
  phi.Redim(1,1);
  dphi.Redim(0,1);
  phi(0,0) = 1.;
}

void TPZCompElPoint::Shape(TPZVec<REAL> &,TPZFMatrix &phi,TPZFMatrix &dphi) {
  phi.Redim(1,1);
  dphi.Redim(0,1);
  phi(0,0) = 1.;
}
//Cedric 16/03/99
void TPZCompElPoint::SetIntegrationRule(int /*order*/) {
}

void TPZCompElPoint::CreateGraphicalElement(TPZGraphMesh &/*grmesh*/, int /*dimension*/) {
}

