// -*- c++ -*-
//METHODS DEFINITION FOR CLASS ELEMGC3D
#include "pzelgc3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgq2d.h"
#include "pzelcc3d.h"
#include "pzshapecube.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZRefCube.h"
#include "TPZGeoCube.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElC3d *gel,TPZCompMesh &mesh,int &index) {
  return new TPZCompElC3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElC3d::fp)(TPZGeoElC3d *,TPZCompMesh &,int &) = CreateEl;

TPZGeoElC3d::TPZGeoElC3d(int id,TPZVec<int> &nodeindexes,int matid,TPZGeoMesh &mesh):
  TPZGeoEl(id,matid,mesh) {
  
  int i,nnod = nodeindexes.NElements();
  if(nnod!=8) {
    PZError << "TPZGeoElC3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<8;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElC3d::TPZGeoElC3d(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=8) {
    PZError << "TPZGeoElC3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<8;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElC3d::~TPZGeoElC3d() {}

void TPZGeoElC3d::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {

	TPZGeoCube::Shape(pt,phi,dphi);
}

TPZGeoElC3d *TPZGeoElC3d::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElC3d(nodeindexes,matind,mesh);
}

int TPZGeoElC3d::NNodes() {
	return TPZGeoCube::NNodes;
}

int TPZGeoElC3d::NodeIndex(int node) {
  if(node<0 || node>7) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElC3d::NSideNodes(int side) {
	return TPZShapeCube::NSideNodes(side);
}

int TPZGeoElC3d::SideNodeIndex(int side,int node) {
	int loc = TPZShapeCube::SideNodeLocId(side,node);
	if(loc<0) return -1;
	return fNodeIndexes[loc];
}

int TPZGeoElC3d::SideNodeLocIndex(int side,int node) {
	return TPZShapeCube::SideNodeLocId(side,node);
}

void TPZGeoElC3d::MidSideNodeIndex(int side,int &index) {
	TPZRefCube::MidSideNodeIndex(this,side,index);
}

void TPZGeoElC3d::NewMidSideNode(int side,int &index) {
	TPZRefCube::NewMidSideNode(this,side,index);
}

int TPZGeoElC3d::SideDimension(int side) {
	return TPZShapeCube::SideDimension(side);
}

/*
TPZGeoElSide TPZGeoElC3d::HigherDimensionSides(int side,int targetdimension) {
  //targetdimension deve ser 1 , 2 ou 3
  //se side =  0 a  7 targetdimension deve ser 1
  //se side =  8 a 19 targetdimension deve ser 2
  //se side = 26      targetdimension deve ser 3
  if( (side<0 || side>25) || (targetdimension < 1 || targetdimension > 3) ) {
    PZError << "TPZGeoElC3d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;//Cedric 26/05/99
  //side = 0 a 25
  switch(targetdimension) {//=1,2
  case 1:
    if(this == father->SubElement(0)) {
      if(side==1) return TPZGeoElSide(this,8);
      if(side==3) return TPZGeoElSide(this,11);
      if(side==4) return TPZGeoElSide(this,12);
    } else if(this == father->SubElement(1)) {
      if(side==0) return TPZGeoElSide(this,8);
      if(side==2) return TPZGeoElSide(this,9);
      if(side==5) return TPZGeoElSide(this,13);
    } else if(this == father->SubElement(2)) {
      if(side==1) return TPZGeoElSide(this,9);
      if(side==3) return TPZGeoElSide(this,10);
      if(side==6) return TPZGeoElSide(this,14);
    } else if(this == father->SubElement(3)) {
      if(side==0) return TPZGeoElSide(this,11);
      if(side==2) return TPZGeoElSide(this,10);
      if(side==7) return TPZGeoElSide(this,15);
    } else if(this == father->SubElement(4)) {
      if(side==0) return TPZGeoElSide(this,12);
      if(side==5) return TPZGeoElSide(this,16);
      if(side==7) return TPZGeoElSide(this,19);
    } else if(this == father->SubElement(5)) {
      if(side==1) return TPZGeoElSide(this,13);
      if(side==4) return TPZGeoElSide(this,16);
      if(side==6) return TPZGeoElSide(this,17);
    } else if(this == father->SubElement(6)) {
      if(side==2) return TPZGeoElSide(this,14);
      if(side==5) return TPZGeoElSide(this,17);
      if(side==7) return TPZGeoElSide(this,18);
    } else if(this == father->SubElement(7)) {
      if(side==3) return TPZGeoElSide(this,15);
      if(side==4) return TPZGeoElSide(this,19);
      if(side==6) return TPZGeoElSide(this,18);
    }
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  case 2:
    if(this == father->SubElement(0)) {

      bestface = BestDimensionSideOfTwoFaces(20,21);
      if((side== 1 || side== 8) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(20,24);
      if((side== 3 || side==11) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(21,24);
      if((side== 4 || side==12) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 2) return TPZGeoElSide(this,20);
      if(side== 9) return TPZGeoElSide(this,20);
      if(side==10) return TPZGeoElSide(this,20);
      if(side== 5) return TPZGeoElSide(this,21);
      if(side==13) return TPZGeoElSide(this,21);
      if(side==16) return TPZGeoElSide(this,21);
      if(side== 7) return TPZGeoElSide(this,24);
      if(side==15) return TPZGeoElSide(this,24);
      if(side==19) return TPZGeoElSide(this,24);
    } else if(this == father->SubElement(1)) {

      bestface = BestDimensionSideOfTwoFaces(20,21);
      if((side== 0 || side== 8) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(20,22);
      if((side== 2 || side== 9) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(21,22);
      if((side== 5 || side==13) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 3) return TPZGeoElSide(this,20);
      if(side==10) return TPZGeoElSide(this,20);
      if(side==11) return TPZGeoElSide(this,20);
      if(side== 4) return TPZGeoElSide(this,21);
      if(side==12) return TPZGeoElSide(this,21);
      if(side==16) return TPZGeoElSide(this,21);
      if(side== 6) return TPZGeoElSide(this,22);
      if(side==14) return TPZGeoElSide(this,22);
      if(side==17) return TPZGeoElSide(this,22);
    } else if(this == father->SubElement(2)) {

      bestface = BestDimensionSideOfTwoFaces(20,22);
      if((side== 1 || side== 9) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(20,23);
      if((side== 3 || side==10) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(22,23);
      if((side== 6 || side==14) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 0) return TPZGeoElSide(this,20);
      if(side== 8) return TPZGeoElSide(this,20);
      if(side==11) return TPZGeoElSide(this,20);
      if(side== 5) return TPZGeoElSide(this,22);
      if(side==13) return TPZGeoElSide(this,22);
      if(side==17) return TPZGeoElSide(this,22);
      if(side== 7) return TPZGeoElSide(this,23);
      if(side==15) return TPZGeoElSide(this,23);
      if(side==18) return TPZGeoElSide(this,23);
    } else if(this == father->SubElement(3)) {

      bestface = BestDimensionSideOfTwoFaces(20,23);
      if((side== 2 || side==10) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(20,24);
      if((side== 0 || side==11) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(23,24);
      if((side== 7 || side==15) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 1) return TPZGeoElSide(this,20);
      if(side== 8) return TPZGeoElSide(this,20);
      if(side== 9) return TPZGeoElSide(this,20);
      if(side== 6) return TPZGeoElSide(this,23);
      if(side==14) return TPZGeoElSide(this,23);
      if(side==18) return TPZGeoElSide(this,23);
      if(side== 4) return TPZGeoElSide(this,24);
      if(side==12) return TPZGeoElSide(this,24);
      if(side==19) return TPZGeoElSide(this,24);
    } else if(this == father->SubElement(4)) {

      bestface = BestDimensionSideOfTwoFaces(21,24);
      if((side== 0 || side==12) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(21,25);
      if((side== 5 || side==16) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(24,25);
      if((side== 7 || side==19) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 1) return TPZGeoElSide(this,21);
      if(side== 8) return TPZGeoElSide(this,21);
      if(side==13) return TPZGeoElSide(this,21);
      if(side== 3) return TPZGeoElSide(this,24);
      if(side==11) return TPZGeoElSide(this,24);
      if(side==15) return TPZGeoElSide(this,24);
      if(side== 6) return TPZGeoElSide(this,25);
      if(side==17) return TPZGeoElSide(this,25);
      if(side==18) return TPZGeoElSide(this,25);
    } else if(this == father->SubElement(5)) {

      bestface = BestDimensionSideOfTwoFaces(21,22);
      if((side== 1 || side==13) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(21,25);
      if((side== 4 || side==16) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(22,25);
      if((side== 6 || side==17) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 0) return TPZGeoElSide(this,21);
      if(side== 8) return TPZGeoElSide(this,21);
      if(side==12) return TPZGeoElSide(this,21);
      if(side== 2) return TPZGeoElSide(this,22);
      if(side== 9) return TPZGeoElSide(this,22);
      if(side==14) return TPZGeoElSide(this,22);
      if(side== 7) return TPZGeoElSide(this,25);
      if(side==18) return TPZGeoElSide(this,25);
      if(side==19) return TPZGeoElSide(this,25);
    } else if(this == father->SubElement(6)) {

      bestface = BestDimensionSideOfTwoFaces(22,23);
      if((side== 2 || side==14) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(22,25);
      if((side== 5 || side==17) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(23,25);
      if((side== 7 || side==18) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 1) return TPZGeoElSide(this,22);
      if(side== 9) return TPZGeoElSide(this,22);
      if(side==13) return TPZGeoElSide(this,22);
      if(side== 3) return TPZGeoElSide(this,23);
      if(side==10) return TPZGeoElSide(this,23);
      if(side==15) return TPZGeoElSide(this,23);
      if(side== 4) return TPZGeoElSide(this,25);
      if(side==16) return TPZGeoElSide(this,25);
      if(side==19) return TPZGeoElSide(this,25);
    } else if(this == father->SubElement(7)) {

      bestface = BestDimensionSideOfTwoFaces(23,24);
      if((side== 3 || side==15) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(23,25);
      if((side== 6 || side==18) && bestface) return TPZGeoElSide(this,bestface);
      bestface = BestDimensionSideOfTwoFaces(24,25);
      if((side== 4 || side==19) && bestface) return TPZGeoElSide(this,bestface);
      if(side== 0) return TPZGeoElSide(this,24);
      if(side==11) return TPZGeoElSide(this,24);
      if(side==12) return TPZGeoElSide(this,24);
      if(side== 2) return TPZGeoElSide(this,23);
      if(side==10) return TPZGeoElSide(this,23);
      if(side==14) return TPZGeoElSide(this,23);
      if(side== 5) return TPZGeoElSide(this,25);
      if(side==16) return TPZGeoElSide(this,25);
      if(side==17) return TPZGeoElSide(this,25);
    }
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  case 3:
    return TPZGeoElSide(this,26);//0<=side<=25
  default:
    cout << "TPZGeoElC3d::Error HigherDimensionSides: targetdimension: " << targetdimension << "\tside " << side << endl;
  }//switch
  return TPZGeoElSide();
}
*/

void TPZGeoElC3d::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides) {
	TPZStack<int> high;
	TPZShapeCube::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}

}

/*
int TPZGeoElC3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

  TPZGeoElSide grandfather = Father(face1);
  int levnum1=0,levnum2=0;
  while(grandfather.Element()) {
    levnum1++;
    grandfather = grandfather.Father();//face1
  }
  grandfather = Father(face2);
  while(grandfather.Element()) {
    levnum2++;
    grandfather = grandfather.Father();//face2
  }
  if(levnum1 >= levnum2) return face1; // "=" Cesar 31/07/02
  if(levnum2 > levnum1) return face2;
  return 0;
}
*/

void TPZGeoElC3d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapeCube::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapeCube::SideConnectLocId(side,is));

}

void TPZGeoElC3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
	TPZFMatrix xco(3,TPZGeoCube::NNodes);
	int i,j;
	for(i=0; i<TPZShapeCube::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoCube::Jacobian(xco,param,jacobian,axes,detjac,jacinv);
}

void TPZGeoElC3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
	TPZFMatrix xco(3,TPZGeoCube::NNodes);
	int i,j;
	for(i=0; i<TPZShapeCube::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoCube::X(xco,loc,result);
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/

/*
void TPZGeoElC3d::NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,
			       TPZFMatrix &axes,TPZFMatrix &jac1d) {

#ifdef DEBUG
  if (nnodes != 8) {
    PZError << "TPZGeoElC3d.NormalVector, only implemented for"
      " 8 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1. || param[2] < -1. || param[2] > 1.) {
    PZError << "TPZGeoElC3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "TPZGeoElC3d::NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side > 5) {//6 faces
    PZError << "TPZGeoElC3d.jacobian invalid side : "
      " side = " << side << "\n";
    return;
  }
#endif

  REAL spacephi[12],spacedphi[30];
  TPZFMatrix phi(8,1,spacephi,12);
  TPZFMatrix dphi(3,8,spacedphi,30);
  Shape(param,phi,dphi);
  TPZGeoNode *np;
  TPZVec<REAL> n(3,0.);
  int i,j,ider;
  if(side==0 || side==5) ider = 2;
  if(side==1 || side==3) ider = 1;
  if(side==2 || side==4) ider = 0;
  for(i=0;i<4;i++) {
    np = NodePtr(TPZCompElC3d::FaceNodes[side][i]);
    for(j=0;j<3;j++)
      n[j] += np->Coord(j)*dphi(ider,i);
  }

  jac1d(0,0) = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

  for(i=0;i<3;i++) normal[i] = n[i]/jac1d(0,0);

  switch(side) {
  case 0: for(i=0;i<3;i++) normal[i] *= -1.;
    break;
  case 1: for(i=0;i<3;i++) normal[i] *= -1.;
    break;
  case 4: for(i=0;i<3;i++) normal[i] *= -1.;
  default:
    cout << "TPZGeoElC3d::X undefined switch case... side = " << side << endl;
  }
  axes.Zero();
  axes(0,0) = 1.0;
  axes(1,1) = 1.0;
  axes(2,2) = 1.0;
  return;
}
*/
/** TO SUBDIVISION
********************************************************************************
Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElC3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
	TPZRefCube::Divide(this,SubElVec);
}


int TPZGeoElC3d::NSubElements() {
	return TPZRefCube::NSubEl;
}

//int TPZGeoElC3d::NSideSubElements(int side) {
//	return TPZRefCube::NSideSubElements(side);
//}

/*
TPZGeoElSide TPZGeoElC3d::SideSubElement(int side,int position) {
  if (position<0 || position>8 || side <0 ||side>26) {
    PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side << endl;
    return TPZGeoElSide();
  }                              //fSubEl[is]
  if(side==26) return TPZGeoElSide(SubElement(position),26);//centro
  if(side<8) {//cantos
    if(position!=0) {
      PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
    } else {
      return TPZGeoElSide(SubElement(side),side);
    }
  }
  if(side>7 && side<20) {//lados
    if(position!=0 && position!=1) {
      PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
    } else {
      int s = side-8;
      return TPZGeoElSide(SubElement(TPZCompElC3d::SideNodes[s][position]),side);
    }
  }
  if(side>19) {//faces
    if(position<0 || position>3) {//position!=0 && position!=1 && position!=2 && position!=3
      PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
    } else {
      int s = side-20;
      return TPZGeoElSide(SubElement(TPZCompElC3d::FaceNodes[s][position]),side);
    }
  }
  return TPZGeoElSide();
}
*/
/*
void TPZGeoElC3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
  if(!fSubEl[0]) {
    sub.Resize(0);
    return;
  }
  if(side < 0 || side > 26) {
    PZError << "TPZGeoElC3d::SideSubElements called for side " << side << endl;
    return;
  }
  if(side==26) {
    sub.Resize(8);
    for(int i=0;i<8;i++) sub[i] = fSubEl[i];
    return;
  }
  if(side<8) {
    sub.Resize(1);
    sub[0]=fSubEl[side];
    return;
  }
  if(side>7 && side<20) {//lados
    int s = side-8;
    sub.Resize(2);
    sub[0] = fSubEl[TPZCompElC3d::SideNodes[s][0]];
    sub[1] = fSubEl[TPZCompElC3d::SideNodes[s][1]];
    return;
  }
  if(side>19) {//faces
    int s = side-20;
    sub.Resize(4);
    sub[0] = fSubEl[TPZCompElC3d::FaceSons[s][0]];
    sub[1] = fSubEl[TPZCompElC3d::FaceSons[s][1]];
    sub[2] = fSubEl[TPZCompElC3d::FaceSons[s][2]];
    sub[3] = fSubEl[TPZCompElC3d::FaceSons[s][3]];
    return;
  }
}
*/
/*
TPZGeoElSide TPZGeoElC3d::Father(int side) {

  if(!fFather) return TPZGeoElSide();
  int whichsub = -1;
  int i;
  for(i=0; i<8; i++) if(fFather->SubElement(i) == this) whichsub = i;
  if(whichsub == -1) {//equivale a is = 4  ou is > 3
    PZError << "TPZGeoElC3d::Father. fFather isn't father of this element.\n";
    return TPZGeoElSide();
  }
  //agora o atual elemento é o filho numero whichsub < 8
  if(side == whichsub || side==26) return TPZGeoElSide(fFather,side);//cantos
  //lados
  if(whichsub == 0 && (side== 8 || side==11 || side==12)) return TPZGeoElSide(fFather,side);
  if(whichsub == 1 && (side== 8 || side== 9 || side==13)) return TPZGeoElSide(fFather,side);
  if(whichsub == 2 && (side== 9 || side==10 || side==14)) return TPZGeoElSide(fFather,side);
  if(whichsub == 3 && (side==10 || side==11 || side==15)) return TPZGeoElSide(fFather,side);
  if(whichsub == 4 && (side==12 || side==16 || side==19)) return TPZGeoElSide(fFather,side);
  if(whichsub == 5 && (side==13 || side==17 || side==16)) return TPZGeoElSide(fFather,side);
  if(whichsub == 6 && (side==14 || side==17 || side==18)) return TPZGeoElSide(fFather,side);
  if(whichsub == 7 && (side==15 || side==18 || side==19)) return TPZGeoElSide(fFather,side);
  //faces
  if(whichsub == 0 && (side==20 || side==21 || side==24)) return TPZGeoElSide(fFather,side);
  if(whichsub == 1 && (side==20 || side==21 || side==22)) return TPZGeoElSide(fFather,side);
  if(whichsub == 2 && (side==20 || side==22 || side==23)) return TPZGeoElSide(fFather,side);
  if(whichsub == 3 && (side==20 || side==23 || side==24)) return TPZGeoElSide(fFather,side);
  if(whichsub == 4 && (side==21 || side==24 || side==25)) return TPZGeoElSide(fFather,side);
  if(whichsub == 5 && (side==21 || side==22 || side==25)) return TPZGeoElSide(fFather,side);
  if(whichsub == 6 && (side==22 || side==23 || side==25)) return TPZGeoElSide(fFather,side);
  if(whichsub == 7 && (side==23 || side==24 || side==25)) return TPZGeoElSide(fFather,side);
  //outro caso
  //  cout <<"TPZGeoElC3d::Father:undefined father for side " << side << endl;
  return TPZGeoElSide();
}
*/
/*
void TPZGeoElC3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

  int nsub = NSideSubElements(side);
  if(!nsub) return;
  sub.Resize(nsub);
  int i,j,k;
  if(nsub==1) {//side = 0 a 7
    if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      PZError << "TPZGeoElC3d::GetSubElement subelement does not contain refnode" << endl;
      return;
    }
    sub[0]=TPZGeoElSide(fSubEl[side],side);
    return;
  }
  //int isub=0;
  for(i=0;i<nsub;i++) {
    TPZGeoElSide sidesub = SideSubElement(side,i);
    TPZGeoEl *subel = sidesub.Element();
    for(k = 0; k < refnode.NElements(); k++) {
      for(j=0;j<8;j++) {
	if(subel->NodeIndex(j)==refnode[k]) {
	  sub[k] = SideSubElement(side,i);
	}
      }
    }
  }
  return;
}
*/
/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
/*
void TPZGeoElC3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform &t) {
  if(!fFather) return;
  int whichsub = -1;
  int i;
  for(i=0; i<8; i++) if(fFather->SubElement(i) == this) whichsub = i;
  if(whichsub == -1) return;
  int dim = SideDimension(side);
  TPZTransform tloc(dim);
  REAL store[12];
  TPZFMatrix mult(dim,dim,store,9);
  TPZFMatrix sum(dim,1,store+9,3);
  mult.Zero();
  sum.Zero();

  TPZGeoEl *locfather;
  if(father == this) {
    mult(0,0) = mult(1,1) = mult(2,2) = 1.;//pai para pai = identidade
    locfather = fFather;
  } else {
    locfather = fFather;//pai do lemento atual
  }
  if(!locfather) {
    cout << "TPZGeoElC3d::BuildTransform could not identify the father element\n";
    return;
  }

  if(side == 26) {//pai para filho
    mult(0,0) = 0.5;
    mult(1,1) = 0.5;
    mult(2,2) = 0.5;
    sum(0,0) = -0.5;
    sum(1,0) = -0.5;
    sum(2,0) = -0.5;
    switch(whichsub) {//o atual é o filho numero whichsub
    case 0:
      break;          //x<0 y<0 z<0
    case 1:
      sum(0,0) *= -1.;//x>0 y<0 z<0
      break;
    case 2:
      sum(0,0) *= -1.;//x>0 y>0 z<0
      sum(1,0) *= -1.;
      break;
    case 3:
      sum(1,0) *= -1.;//x<0 y>0 z<0
      break;
    case 4:
      sum(2,0) *= -1.;//x<0 y<0 z>0
      break;
    case 5:
      sum(0,0) *= -1.;//x>0 y<0 z>0
      sum(2,0) *= -1.;
      break;
    case 6:
      sum(0,0) *= -1.;//x>0 y>0 z>0
      sum(1,0) *= -1.;
      sum(2,0) *= -1.;
      break;
    case 7:
      sum(1,0) *= -1.;//x<0 y>0 z>0
      sum(2,0) *= -1.;
      break;
    default:
      cout << "TPZGeoElC3d::BuildTransform could not identify which sub for side 26 - whichsub = " << whichsub << endl;
    }
  } else
    if(side>19) {//face do pai para face do filho
      whichsub = -1;
      int s = side-20;
      for(i=0; i<4; i++) if(fFather->SubElement(TPZCompElC3d::FaceSons[s][i]) == this) whichsub = i;
      if(whichsub == -1) return;
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      sum(0,0) = -0.5;
      sum(1,0) = -0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
      case 0:
	break;
      case 1:
	sum(0,0) *= -1.;
	break;
      case 2:
	sum(0,0) *= -1.;
	sum(1,0) *= -1.;
	break;
      case 3:
	sum(1,0) *= -1.;
	break;
      default:
	cout << "TPZGeoElC3d::BuildTransform could not identify which sub for side " <<  side << " - whichsub = " << whichsub << endl;
      }
    } else
      if(side>7) {//e side < 20 : arestas
	whichsub = -1;
	int s = side-8;//0,1,2,3,4,5,6,7,8,9,10,11
	for(i=0; i<2; i++)
	  if(fFather->SubElement(TPZCompElC3d::SideNodes[s][i]) == this) whichsub = i;
	if(whichsub == -1) return;
   	mult(0,0) = 0.5;
	if(whichsub==0) sum(0,0) = -0.5;
	if(whichsub==1) sum(0,0) =  0.5;
      }
  tloc.SetMatrix(mult,sum);
  t = tloc.Multiply(t);
  if(locfather != father) locfather->BuildTransform(side,father,t);
}
*/

TPZTransform TPZGeoElC3d::SideToSideTransform(int sidefrom,int sideto) {
	return TPZShapeCube::SideToSideTransform(sidefrom,sideto);
}

TPZGeoEl *TPZGeoElC3d::CreateBCGeoEl(int side,int bc) {
	TPZGeoEl *gel = TPZGeoCube::CreateBCGeoEl(this,side,bc);
	return gel;
}

/*void TPZGeoElC3d::NodeFaceIds(TPZVec<int> &ids,int face) {

  ids.Resize(4,-1);
  switch(face) {
  case 20:
    ids[0] = NodeIndex(0);
    ids[1] = NodeIndex(1);
    ids[2] = NodeIndex(2);
    ids[3] = NodeIndex(3);
    return;
  case 21:
    ids[0] = NodeIndex(0);
    ids[1] = NodeIndex(1);
    ids[2] = NodeIndex(5);
    ids[3] = NodeIndex(4);
    return;
  case 22:
    ids[0] = NodeIndex(1);
    ids[1] = NodeIndex(2);
    ids[2] = NodeIndex(6);
    ids[3] = NodeIndex(5);
    return;
  case 23:
    ids[0] = NodeIndex(3);
    ids[1] = NodeIndex(2);
    ids[2] = NodeIndex(6);
    ids[3] = NodeIndex(7);
    return;
  case 24:
    ids[0] = NodeIndex(0);
    ids[1] = NodeIndex(3);
    ids[2] = NodeIndex(7);
    ids[3] = NodeIndex(4);
    return;
  case 25:
    ids[0] = NodeIndex(4);
    ids[1] = NodeIndex(5);
    ids[2] = NodeIndex(6);
    ids[3] = NodeIndex(7);
    return;
  default :
    cout << "TPZCompElC3d::NodeFaceIds bad side , side = " << face << endl;
  }
}
*/

TPZGeoElSide TPZGeoElC3d::Father2(int side){
	if(!fFather) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0) return TPZGeoElSide();
	return TPZGeoElSide(fFather,TPZRefCube::FatherSide(side,son));
}


int TPZGeoElC3d::FatherSide(int side, int son){
	return TPZRefCube::FatherSide(side,son);
}

void TPZGeoElC3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

	TPZRefCube::GetSubElements(this,side,subel);
}

int TPZGeoElC3d::NSideSubElements2(int side) {
	return TPZRefCube::NSideSubElements(side);
}


TPZTransform TPZGeoElC3d::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t){
	
	if(side<0 || side>(TPZShapeCube::NSides-1) || !fFather){
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

static REAL MidSideNode[27][3] = {
  /*00*/{-1.,-1.,-1.},/*01*/{1.,-1.,-1.},/*02*/{1.,1.,-1.},/*03*/{-1.,1.,-1.},
  /*04*/{-1.,-1., 1.},/*05*/{1.,-1., 1.},/*06*/{1.,1., 1.},/*07*/{-1.,1., 1.},
  /*08*/{ 0.,-1.,-1.},/*09*/{1., 0.,-1.},/*10*/{0.,1.,-1.},/*11*/{-1.,0.,-1.},
  /*12*/{-1.,-1., 0.},/*13*/{1.,-1., 0.},/*14*/{1.,1., 0.},/*15*/{-1.,1., 0.},
  /*16*/{ 0.,-1., 1.},/*17*/{1., 0., 1.},/*18*/{0.,1., 1.},/*19*/{-1.,0., 1.},
  /*20*/{ 0., 0.,-1.},/*21*/{0.,-1., 0.},/*22*/{1.,0., 0.},/*23*/{ 0.,1., 0.},
  /*24*/{-1., 0., 0.},/*25*/{0., 0., 1.},/*26*/{0.,0., 0.} };

int TPZGeoElC3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
  //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<8;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<27;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      TPZTransform telsd = TPZShapeCube::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> son side 
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform tson(son->SideDimension(sd));
      TPZTransform t = son->BuildTransform2(sd,gel,tson);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = TPZShapeCube::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(26).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao errada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(26).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(26).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}

void TPZGeoElC3d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >7){
    PZError << "TPZGeoElC3d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}
   

TPZIntPoints * TPZGeoElC3d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoCube::CreateSideIntegrationRule(side,order);
}

TPZTransform TPZGeoElC3d::GetTransform(int side,int son) {
	return TPZRefCube::GetTransform(side,son);
}

void TPZGeoElC3d::CenterPoint(int side, TPZVec<REAL> &masscent){
  TPZShapeCube::CenterPoint(side,masscent);
}

REAL TPZGeoElC3d::RefElVolume(){
  return TPZShapeCube::RefElVolume();
}
