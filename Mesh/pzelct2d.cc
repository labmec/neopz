//$Id: pzelct2d.cc,v 1.9 2003-11-06 19:15:19 cesar Exp $

// -*- c++ -*-
#include "pzelct2d.h"
#include "pzelcq2d.h"
#include "pzerror.h"
//#include "pzelgt2d.h"
#include "pzgeoel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzconnect.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzmat2dlin.h"
#include "pzcmesh.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
//#include "pzgraphmesh.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "pzgraphel.h"
//#include "pzgraphnode.h"


TPZCompElT2d::TPZCompElT2d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZInterpolatedElement(mesh,ref,index), fIntRule(1) {

  int i;
  for(i=0; i<7; i++) fConnectIndexes[i]=-1;
  for(i=0;i<4;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
  // Jorge 19/5/99
  /** Compatibilizing the continuity over the neighboard elements */
  for(i=0;i<7;i++) {
    TPZStack<TPZCompElSide> elvec;
    TPZCompElSide thisside(this,i);
    thisside.EqualLevelElementList(elvec,0,0);
    /** Asking if exist neighbour : If this neighboard is continuous all neighboards are continuous */
    if(elvec.NElements()) {
      ((TPZInterpolatedElement *)elvec[0].Element())->MakeConnectContinuous(elvec[0].Side());
    }
  }

//  RemoveSideRestraintsII(EInsert);
  ref->SetReference(this);
  for(i=0;i<7;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(i=3;i<7;i++)
    IdentifySideOrder(i);
  TPZVec<int> order(2,2*fSideOrder[3]);
  fIntRule.SetOrder(order);
}

TPZCompElT2d::TPZCompElT2d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index,int /*noconnects*/) :
		TPZInterpolatedElement(mesh,ref,index), fIntRule(1) {
  int i;
  for(i=0; i<7; i++) fConnectIndexes[i]=-1;
  for(i=0;i<4;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
//  RemoveSideRestraintsII(EInsert);
}

TPZCompElT2d::TPZCompElT2d(TPZCompMesh &mesh, const TPZCompElT2d &copy) :
		TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
  int i;
  for(i=0; i<7; i++) fConnectIndexes[i]=copy.fConnectIndexes[i];
  for(i=0;i<4;i++) {
    fSideOrder[i] = copy.fSideOrder[i];
    fPreferredSideOrder[i] = copy.fPreferredSideOrder[i];
  }
}

void TPZCompElT2d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<7) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElT2d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElT2d::ConnectIndex(int i) {
	if(i<0 || i>6) {
   	PZError << "TPZCompElT2d::ConnectIndex. Bad parameter i.\n";
      return -1;
   }
   return fConnectIndexes[i];
}

int TPZCompElT2d::NConnectShapeF(int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return 1;
  case 3:
  case 4:
  case 5:
    return fSideOrder[side-3]-1;
  case 6:
    return (fSideOrder[3]-2) < 0 ? 0 : ((fSideOrder[3]-2)*(fSideOrder[3]-1))/2;
  default:
    PZError << "TPZCompElT2d::NConnectShapeF, bad parameter iconnect " << side << endl;
    return 0;
  }
}

int TPZCompElT2d::NSideConnects(int side) {
	return TPZShapeTriang::NSideConnects(side);
}

/**It do not verify the values of the c*/
int TPZCompElT2d::SideConnectLocId(int c, int side) {
	return TPZShapeTriang::SideConnectLocId(side,c);
}

void TPZCompElT2d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=4)
    PZError << "TPZCompElT2d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<4;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElT2d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(4);
  for(int i=0;i<4;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElT2d::CreateSideIntegrationRule(int side) {  // or TPZInt...(2*ord-1) ???
   if(side<3 || side>6) return new TPZInt1d(0);
   if(side<6) return new TPZInt1d(2*SideOrder(side));
   return new TPZIntTriang(2*fSideOrder[3]);
}

int TPZCompElT2d::PreferredSideOrder(int side) {
  if(side<3) return 0;
  if(side<7) {
	  int order = fPreferredSideOrder[side-3];
	  return AdjustPreferredSideOrder(side,order);
  }
  PZError << "TPZCompElT2d::PreferredSideOrder called for side = " << side << "\n";
  return 0;
}

void TPZCompElT2d::SetPreferredSideOrder(int order) {
  int side;
  for(side=3; side<7; side++) fPreferredSideOrder[side-3] = order;
}


void TPZCompElT2d::SetSideOrder(int side, int order) {
/*   if(fConnectIndexes[side] ==-1) { */
/*     cout << "TPZCompElT2d::SetSideOrder called for uninitialized connect\n"; */
/*     return; */
/*   } */
  if(side<0 || side>6 || order <1){
    PZError << "TPZCompElT2d::SetSideOrder. Bad parameter side.\n";
    return;
  }
  if(side>2) fSideOrder[side-3] = order;
  if(fConnectIndexes[side] !=-1) { 
    TPZConnect &c = Connect(side);
    c.SetOrder(order);
    int seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial *mat = Material();
    if(mat) nvar = mat->NStateVariables();
    Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
    if (side == NConnects()-1){
      TPZVec<int> ord(2,2*order+2);
      fIntRule.SetOrder(ord);
    }
  }
}

int TPZCompElT2d::SideOrder(int side) {
  if(side<3 || side>6) return 0;
  return fSideOrder[side-3];
}

void TPZCompElT2d::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {
  point[0] = point[1] = 0.;
  switch(side) {
  case 0:
    return;
  case 1:
    point[0]=1.;
    return;
  case 2:
    point[1]=1.;
    return;
  case 3:
    point[0] = (1.+par[0])*.5;
	 return;
  case 4:
    point[0] = (1.-par[0])*.5;
	 point[1] = (1.+par[0])*.5;
    return;
  case 5:
    point[1] = (1.-par[0])*.5;
    return;
  case 6:
    point[0]=par[0];
    point[1]=par[1];
    return;
  default:
    PZError << "TPZCompElT2d::SideParameterToElement. Bad paramenter side.\n";
  }
}

void TPZCompElT2d::ElementToSideParameter(int side,TPZVec<REAL> &point,TPZVec<REAL> &par) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return;
  case 3:
    par[0] = 2.*point[0]-1.;
    return;
  case 4:
    par[0] = point[1]-point[0];
    return;
  case 5:
    par[0] = 1.- 2*point[1];
    return;
  case 6:
    par[0]=point[0];
    par[1]=point[1];
    return;
  default:
    PZError << "TPZCompElT2d::ElementToSideParameter. Bad side = " << side << ".\n";
  }
}

void TPZCompElT2d::CornerShape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
  TPZShapeTriang::ShapeCorner(pt,phi,dphi);
}

void TPZCompElT2d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {  //    id ???
  if(side<0 || side>6) PZError << "TPZCompElT2d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==6) Shape(point,phi,dphi);
  else if(side<3) {
  	phi(0,0) = 1.;
  } else {
    TPZVec<int> id(2);
    TPZManVector<int,1> sideord(1,SideOrder(side));
    id[0] = fReference->NodePtr(side-3)->Id();
    id[1] = fReference->NodePtr((side-2)%3)->Id();
    TPZShapeLinear::Shape(point,id,sideord,phi,dphi);
  }
}

void TPZCompElT2d::Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi) {
  TPZVec<int> id(3);

  int i;
  for(i=0;i<3;i++) id[i] = fReference->NodePtr(i)->Id();//>NodeIndex(i);
  TPZManVector<int> ord(4);
  for(i=0; i<4; i++) ord[i] = fSideOrder[i];
  TPZShapeTriang::Shape(x,id,ord,phi,dphi);
}
//Cedric 16/03/99
void TPZCompElT2d::SetIntegrationRule(int order) {
	TPZIntTriang inttriang(order);
   SetIntegrationRule(inttriang);
}

void TPZCompElT2d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
	  if(dimension == 2) new TPZGraphElTd(this,&grmesh);
}

TPZCompMesh *TPZCompElT2d::CreateMesh() {
   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("A simple square");
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0;
   coord[1] = 0;
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 3;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 3;
   coord[1] = 3;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0;
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(3);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
    //elementos geometricos
   int index;
   TPZGeoEl *elg0 = firstmesh->CreateGeoElement(ETriangle,nodeindexes,1,index);
   //   TPZGeoEl *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 0;
   nodeindexes[1] = 2;
   nodeindexes[2] = 3;
   //elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);

   firstmesh->BuildConnectivity();

   TPZVec<TPZGeoEl *> sub;
   elg0->Divide(sub);
   sub[0]->Divide(sub);
   TPZGeoElBC(elg0,4,-1,*firstmesh);
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("A simple computational mesh");

   TPZMat2dLin *mat = new TPZMat2dLin(1);
   TPZFMatrix xk(1,1,1.),xc(1,1,1.),xf(1,1,1.);
   //   xk(0,1) = xk(1,0) = xc(0,1) = xc(1,0) = 0.;
   mat->SetMaterial(xk,xc,xf);
   TPZFMatrix val1(1,1,0.),val2(1,1,0.);
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   secondmesh->InsertMaterialObject(mat);
   secondmesh->InsertMaterialObject(bc);
   secondmesh->AutoBuild();



   return secondmesh;
}
