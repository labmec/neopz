// -*- c++ -*-

//METHODS DEFINITION FOR CLASS ELEM1D

#include "pzelc1d.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzconnect.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzshtmat.h"
#include "pzelmat.h"
#include "pzerror.h"
//#include "pzgraphel.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztempmat.h"
#include "pzmaterial.h"
#include "pzshapelinear.h"
#include <math.h>


TPZCompEl1d::TPZCompEl1d(TPZCompMesh &mesh,TPZGeoEl *ref, int &index)
  : TPZInterpolatedElement(mesh,ref,index), fIntRule(1) {

  fSideOrder = gOrder;
  fPreferredSideOrder = gOrder;
//  RemoveSideRestraintsII(EInsert);
  int i;

  for(i=0; i<3; i++) {
    fConnectIndexes[i] = -1;
  }
  /** Jorge 19/5/99 Compatibilizing the continuity over the neighboard elements */
  for(i=0;i<3;i++) {
    TPZStack<TPZCompElSide> elvec;
    TPZCompElSide thisside(this,i);
    thisside.EqualLevelElementList(elvec,0,0);
    /** Asking if exist neighbour. If this neighboard is continuous all neighboards are continuous */
    if(elvec.NElements()) {
      ((TPZInterpolatedElement *)elvec[0].Element())->MakeConnectContinuous(elvec[0].Side());
    }
  }

  ref->SetReference(this);
  for(i=0;i<3;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }

  IdentifySideOrder(2);

  TPZVec<int> order(1,fSideOrder*2);
  order [0] = 20;
  fIntRule.SetOrder(order);
}

TPZCompEl1d::TPZCompEl1d(TPZCompMesh &mesh,TPZGeoEl *ref, int &index, int)
  : TPZInterpolatedElement(mesh,ref,index), fIntRule(1) {

  int i;
  for(i=0; i<3; i++) fConnectIndexes[i] = -1;
  fSideOrder = gOrder;
  fPreferredSideOrder = gOrder;
//  RemoveSideRestraintsII(EInsert);
}

TPZCompEl1d::~TPZCompEl1d(){
  	if(Reference()) {
      if(Reference()->Reference()) {
         RemoveSideRestraintsII(EDelete);
      }
      Reference()->ResetReference();
   }
}



void TPZCompEl1d::SetConnectIndex(int i,int connectindex) {
  if(i>=0 && i<3) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompEl1d::SetConnectIndex. Bad parameter i = " << i << " .\n";
    PZError.flush();
  }
}

int TPZCompEl1d::ConnectIndex(int i) {
  if(i>=0 && i<3) return fConnectIndexes[i];

  PZError << "TPZCompEl1d::ConnectIndex. Bad parameter i = " << i << " .\n";
  PZError.flush();
  return -1;
}

int TPZCompEl1d::NSideConnects(int i) {
	return TPZShapeLinear::NSideConnects(i);
}

int TPZCompEl1d::NConnectShapeF(int iconnect) {
  if(iconnect == 1 || iconnect == 0) return 1;
  if(iconnect == 2) return fSideOrder-1;
  PZError << "TPZCompEl1d::NConnectShapeF, bad parameter iconnect " << iconnect << endl;
  return 0;
}

int TPZCompEl1d::SideConnectLocId(int c,int side) {
	return TPZShapeLinear::SideConnectLocId(side,c);
}

void TPZCompEl1d::Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi) {
  TPZVec<int> id(2);
  TPZManVector<int,1> sideorder(1,fSideOrder);
  id[0] = fReference->NodePtr(0)->Id();
  id[1] = fReference->NodePtr(1)->Id();
  TPZShapeLinear::Shape(x,id,sideorder,phi,dphi);  //Shape1d(x[0],NumberOfShapeF(),phi,dphi,id);
}            //Porque utilizar os ids dos elementos geometricos id[1] < id[0] ?

void TPZCompEl1d::SetInterpolationOrder(TPZVec<int> &ord) {
  fPreferredSideOrder = ord[0];
}

void TPZCompEl1d::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(1);
	ord[0] = fSideOrder;
}

int TPZCompEl1d::PreferredSideOrder(int side) {
  if(side == 0 || side == 1) return 0;
  if(side == 2) return fPreferredSideOrder;
  return -1;
}

void TPZCompEl1d::SetPreferredSideOrder(int side, int order) {
	if(side < 2) return;
   if(side < 3) fPreferredSideOrder = order;
   return;
}

void TPZCompEl1d::SetSideOrder(int side, int order) {
  //  if(fConnectIndexes[side] ==-1) {
  //cout << "TPZCompEl1d::SetSideOrder called for uninitialized connect\n";
  //  return;
  //}
   if(side == 2) {
	   if(order < 1) return;
      fSideOrder = order;
       if(fConnectIndexes[side] !=-1) {
	 TPZConnect &c = Connect(side);
	 int seqnum = c.SequenceNumber();
	 int nvar = 1;
	 TPZMaterial *mat = Material();
	 if(mat) nvar = mat->NStateVariables();
	 Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
	 TPZVec <int> ord(1,2*order+2);
	 fIntRule.SetOrder(ord);
       }
   } else if(order)
   	PZError << "TPZCompEl1d::SetSideOrder side = " << side << " order = " << order << endl;
}

int TPZCompEl1d::SideOrder(int side) {
  if(side == 0 || side == 1) return 0;
  if(side == 2) return fSideOrder;
  return -1;
}

void TPZCompEl1d::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {
  if(side == 0) point[0] = -1.;
  else if(side == 1) point[0] = 1;
  else {
    point[0] = par[0];
  }
}

void TPZCompEl1d::ElementToSideParameter(int side,TPZVec<REAL> &point,TPZVec<REAL> &par) {
  if(side == 2) par[0] = point[0];
}

void TPZCompEl1d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {
  if(side == 0 || side == 1) {
    phi(0,0) = 1.;
    //    dphi(0,0) = 0.;
  } else Shape(point,phi,dphi);
}

TPZIntPoints *TPZCompEl1d::CreateSideIntegrationRule(int side){
  if(side != 2) return new TPZInt1Point();
  return new TPZInt1d(2*fSideOrder);
}
//Cedric 16/03/99 
void TPZCompEl1d::SetIntegrationRule(int order) {
	TPZInt1d int1d(order);
   SetIntegrationRule(int1d);
}

void TPZCompEl1d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  if(dimension == 1) new TPZGraphEl1dd(this,&grmesh);
}


