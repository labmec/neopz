//$Id: pzelcq2d.cc,v 1.13 2003-12-12 19:59:20 phil Exp $

// -*- c++ -*-
#include "pzelcq2d.h"
#include "pzfmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
//#include "pzelgq2d.h"
#include "pzcmesh.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzgraphel.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzgraphelq2dd.h"

#include "pzmat2dlin.h"
#include "pzelasmat.h"
//#include "../graf/grafgrid.h"
//#include "../graf/grafel.h"
//#include "../graf/grafnode.h"


TPZCompElQ2d::TPZCompElQ2d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZInterpolatedElement(mesh,ref,index), fIntRule(2,2) {
  int i;
  for(i=0;i<5;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
  for(i=0; i<9; i++) fConnectIndexes[i]=-1;
//  RemoveSideRestraintsII(EInsert);
  ref->SetReference(this);
  for(i=0;i<4;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<9;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
	 IdentifySideOrder(i);
  }

  TPZVec<int> order(2,2*fSideOrder[4]+2);
  //  TPZVec<int> order(2,1);
  fIntRule.SetOrder(order);
}

TPZCompElQ2d::TPZCompElQ2d(TPZCompMesh &mesh,TPZGeoEl *ref,int &index,int /*noconnects*/) :
			TPZInterpolatedElement(mesh,ref,index), fIntRule(2,2) {
  int i;
  for(i=0; i<9; i++) fConnectIndexes[i] = -1;
  for(i=0;i<5;i++) {
    fSideOrder[i] = gOrder;
    fPreferredSideOrder[i] = gOrder;
  }
//  RemoveSideRestraintsII(EInsert);
  TPZVec<int> order(2,2*fSideOrder[4]+2);
  //  TPZVec<int> order(2,1);
  fIntRule.SetOrder(order);
}

TPZCompElQ2d::TPZCompElQ2d(TPZCompMesh &mesh, const TPZCompElQ2d &copy) :
			TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
  int i;
  for(i=0; i<9; i++) fConnectIndexes[i] = copy.fConnectIndexes[i];
  for(i=0;i<5;i++) {
    fSideOrder[i] = copy.fSideOrder[i];
    fPreferredSideOrder[i] = copy.fPreferredSideOrder[i];
  }
}

void TPZCompElQ2d::VarRange(int var,double &min,double &max) {
	PZError << "TPZCompElQ2d::VarRange is not defined.\n";
   if(var>-1) max = min = 0.;
}

void TPZCompElQ2d::Load() {
	PZError << "TPZCompElQ2d::Load is called.\n";
}

void TPZCompElQ2d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<9) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElQ2d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElQ2d::ConnectIndex(int i) {
	if(i>8 || i<0) {
   	PZError << "TCompElT2d::ConnectIndex. Bad parameter i.\n";
      return -1;
   }
   return fConnectIndexes[i];
}

int TPZCompElQ2d::NConnectShapeF(int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return 1;
  case 4:
  case 5:
  case 6:
  case 7:
	  return fSideOrder[side-4]-1;
  case 8:
     return (fSideOrder[side-4]-1)*(fSideOrder[side-4]-1);//Cedric Modified 24/02/99
  default:
    PZError << "TPZCompElQ2d::NConnectShapeF, bad parameter iconnect " << side << endl;
    return 0;
  }
}

int TPZCompElQ2d::NSideConnects(int side) {
	return TPZShapeQuad::NSideConnects(side);
}

/**It do not verify the values of the c*/
int TPZCompElQ2d::SideConnectLocId(int c,int side) {
	return TPZShapeQuad::SideConnectLocId(side,c);
}

void TPZCompElQ2d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=5)
    PZError << "TPZCompElT2d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<5;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElQ2d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(5);
  for(int i=0;i<5;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElQ2d::CreateSideIntegrationRule(int side) {  // or TPZInt...(2*ord-1) ???

   if(side<4 || side>8) return new TPZInt1Point();
   if(side<8) return new TPZInt1d(2*SideOrder(side));
   return new TPZIntQuad(2*fSideOrder[4],2*fSideOrder[4]);
}

int TPZCompElQ2d::PreferredSideOrder(int side) {
  if(side<4) return 1;
  if(side<9) {
	  int order =fPreferredSideOrder[side-4];
	  return AdjustPreferredSideOrder(side,order);
  }
   PZError << "TPZCompElQ2d::PreferredSideOrder called for side = " << side << "\n";
  return 0;
}

void TPZCompElQ2d::SetPreferredSideOrder(int order) {
  int side;
  for(side=0; side<5; side++) fPreferredSideOrder[side] = order;
}

void TPZCompElQ2d::SetSideOrder(int side, int order) {
/*   if(fConnectIndexes[side] ==-1) { */
/*     cout << "TPZCompElQ2d::SetSideOrder called for uninitialized connect\n"; */
/*     return; */
/*   } */
  if(side<0 || side>8 || (order < 1 && side > 3)){
    PZError << "TPZCompElQ2d::SetSideOrder. Bad paramenter side.\n" << side << " order = " << order << endl;
    return;
  }
  if(side>3) fSideOrder[side-4] = order;
  if(fConnectIndexes[side] != -1) {
    TPZConnect &c = Connect(side);
    if(c.HasDependency() && order != c.Order()) {
      cout << "TPZCompElQ2d::SetSideOrder fodeu\n";
    }
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

int TPZCompElQ2d::SideOrder(int side) {
  if(side<4 || side>8) return 0;
  return fSideOrder[side-4];
}

void TPZCompElQ2d::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {
  switch(side) {
  case 0:
    point[0]=point[1]=-1.;
    return;
  case 1:
    point[0]= 1.;
    point[1]=-1.;
    return;
  case 2:
    point[0]=point[1]=1.;
    return;
  case 3:
    point[0]=-1.;
    point[1]=1.;
    return;
  case 4:
  case 6:
    point[0] = par[0];
	 point[1] = (side==4) ? -1.0 : 1.0;
    return;
  case 5:
  case 7:
    //point[0] = ( side == 1 ) ? 1.0 : -1.0;
    point[0] = ( side == 5 ) ? 1.0 : -1.0;//Cedric 04/05/99
    point[1] = par[0];
    return;
  case 8:
    point[0]=par[0];
    point[1]=par[1];
    return;
  default:
    PZError << "TPZCompElQ2d::SideParameterToElement. Bad paramenter side.\n";
    return;
  }
}

void TPZCompElQ2d::ElementToSideParameter(int side,TPZVec<REAL> &point,TPZVec<REAL> &par) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return;
  case 4:
    par[0] = point[0];
    return;
  case 6:
    par[0] = -point[0];
    return;
  case 5:
    par[0] = point[1];
    return;
  case 7:
    par[0] = -point[1];
    return;
  case 8:
    par[0]=point[0];
    par[1]=point[1];
    return;
  default:
    PZError << "TCompElQ2d::ElementToSideParameter. Bad paramenter side.\n";
    return;
  }
}

void TPZCompElQ2d::CornerShape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
  TPZShapeQuad::ShapeCornerQuad(pt,phi,dphi);
}

void TPZCompElQ2d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {  //    id ???
  if(side<0 || side>8) PZError << "TPZCompElQ2d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==8) Shape(point,phi,dphi);
  else if(side<4) {
  		phi(0,0)=1.;
  } else {
    TPZVec<int> id(2);
    TPZManVector<int,1> sideord(1,SideOrder(side));
    //    id[0] = fReference->NodeIndex(side-4);
    //    id[1] = fReference->NodeIndex((side-3)%4);
    id[0] =  fReference->NodePtr(side-4)->Id();
    id[1] =  fReference->NodePtr((side-3)%4)->Id();
    TPZShapeLinear::Shape(point,id,sideord,phi,dphi);
  }
}

void TPZCompElQ2d::Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi) {
  TPZVec<int> id(4);//Cedric
  int i;
  for(i=0;i<4;i++) id[i] = fReference->NodePtr(i)->Id();//  fReference->NodeIndex(i);

  TPZManVector<int> ord(5);//Cedric
  for(i=0; i<5; i++) ord[i] = fSideOrder[i];
  TPZShapeQuad::Shape(x,id,ord,phi,dphi);
}

void TPZCompElQ2d::SetIntegrationRule(int order) {
  TPZIntQuad intquad(order,order);
  SetIntegrationRule(intquad);
}

void TPZCompElQ2d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  if(dimension == 2 && Material()->Id() > 0) new TPZGraphElQ2dd(this,&grmesh);
}

void TPZCompElQ2d::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {

  if(var == -1) {
    //    sol[0] = SideOrder(8);
    
    if(fabs(qsi[0]+1.) < 1.e-6) {
      // I am on the left side
      if(fabs(qsi[1]+1.) < 1.e-6) {
	// I am on the bottom side
	sol[0] = (SideOrder(4)+SideOrder(7))/2.;
      } else if(fabs(qsi[1]-1.) < 1.e-6) {
	// I am on the top side
	sol[0] = (SideOrder(6)+SideOrder(7))/2.;
      } else {
	// I am in the middle
	sol[0] = SideOrder(7);
      }
    } else if(fabs(qsi[0]-1.) < 1.e-6) {
      // I am on the right side
      if(fabs(qsi[1]+1.) < 1.e-6) {
	// I am on the bottom side
	sol[0] = (SideOrder(4)+SideOrder(5))/2.;
      } else if(fabs(qsi[1]-1.) < 1.e-6) {
	// I am on the top side
	sol[0] = (SideOrder(5)+SideOrder(6))/2.;
      } else {
	// I am in the middle
	sol[0] = SideOrder(5);
      }
    } else {
      // I am in the middle (horizontally)
      if(fabs(qsi[1]+1.) < 1.e-6) {
	// I am on the bottom side
	sol[0] = SideOrder(4);
      } else if(fabs(qsi[1]-1.) < 1.e-6) {
	// I am on the top side
	sol[0] = SideOrder(6);
      } else {
	// I am in the middle
	sol[0] = SideOrder(8);
      }
    }
  } else {
    TPZInterpolatedElement::Solution(qsi,var,sol);
  }
}

TPZCompMesh *TPZCompElQ2d::CreateMesh() {
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
   TPZVec<int> nodeindexes(4);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   nodeindexes[3] = 3;
    //elementos geometricos
   int index;
   TPZGeoEl *elg0 = firstmesh->CreateGeoElement(EQuadrilateral,nodeindexes,1,index);
   //   TPZGeoEl *elg0 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);

   firstmesh->BuildConnectivity();

   TPZVec<TPZGeoEl *> sub;  elg0->Divide(sub);
   sub[0]->Divide(sub);
   TPZGeoElBC(elg0,4,-1,*firstmesh);
   TPZGeoElBC(elg0,2,-1,*firstmesh);
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("A simple computational mesh");

//    TPZMat2dLin *mat2d = new TPZMat2dLin(1);
//    TPZFMatrix xk(1,1,1.),xc(1,1,1.),xf(1,1,1.);
//    //   xk(0,1) = xk(1,0) = xc(0,1) = xc(1,0) = 0.;
//    mat->SetMaterial(xk,xc,xf);
   TPZElasticityMaterial *mat = new TPZElasticityMaterial(1,2.,0.3,1.,1.);
   TPZFMatrix val1(2,2,0.),val2(2,1,0.);
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   secondmesh->InsertMaterialObject(mat);
   secondmesh->InsertMaterialObject(bc);
   secondmesh->AutoBuild();

   secondmesh->AdjustBoundaryElements();

   return secondmesh;
}

int TPZCompElQ2d::CheckElementConsistency() {

  int n = NCornerConnects();
  int nc = NConnects();
  for(;n<nc; n++) {
    TPZConnect &c = Connect(n);
    int order = c.Order();
    int order2 = fSideOrder[n-NCornerConnects()];
    int ndof = c.NDof(*Mesh());
    int ndof2 = Material()->NStateVariables()*NConnectShapeF(n);
    if(order != order2 || ndof != ndof2) {
      cout << "TPCompElC3d inconsistent datastructure connect order " << order << " side order " << order2 << " blocksize " << ndof <<
	" computed blocksize " << ndof2 << endl;
      Print();
      return 0;
    }
  }
  return 1;
}
