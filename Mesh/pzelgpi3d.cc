//$Id: pzelgpi3d.cc,v 1.4 2003-11-05 16:02:21 tiago Exp $

// -*- c++ -*-
//METHODS DEFINITION FOR CLASS ELEMPI3D
#include "pzelgpi3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelgt3d.h"
#include "pzshapepiram.h"
#include "pzshapetetra.h"
#include "pzelcpi3d.h"
#include "pzelct3d.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElPi3d *gel,TPZCompMesh &mesh,int &index) {
  return new TPZCompElPi3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElPi3d::fp)(TPZGeoElPi3d *,TPZCompMesh &,int &) = CreateEl;

TPZGeoElPi3d::TPZGeoElPi3d(int id,TPZVec<int> &nodeindexes,int matid,TPZGeoMesh &mesh):
  TPZGeoEl(id,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=5) {
    PZError << "TPZGeoElPi3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<5;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<10;i++) fSubEl[i] = 0;
}

TPZGeoElPi3d::TPZGeoElPi3d(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=5) {
    PZError << "TPZGeoElPi3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<5;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<10;i++) fSubEl[i] = 0;
}

TPZGeoElPi3d::~TPZGeoElPi3d() {}

void TPZGeoElPi3d::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
	TPZGeoPyramid::Shape(pt,phi,dphi);
}

TPZGeoElPi3d *TPZGeoElPi3d::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElPi3d(nodeindexes,matind,mesh);
}

int TPZGeoElPi3d::NNodes() {
	return TPZShapePiram::NNodes;
}
int TPZGeoElPi3d::NodeIndex(int node) {
  if(node<0 || node>4) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElPi3d::NSideNodes(int side) {
	return TPZShapePiram::NSideNodes(side);
}

int TPZGeoElPi3d::SideNodeIndex(int side,int node) {
	int loc = TPZShapePiram::SideNodeLocId(side,node);
	if(loc<0) return loc;
	return fNodeIndexes[loc];
}

int TPZGeoElPi3d::SideNodeLocIndex(int side,int node) {
	return TPZShapePiram::SideNodeLocId(side,node);
}

void TPZGeoElPi3d::MidSideNodeIndex(int side,int &index) {
	TPZRefPyramid::MidSideNodeIndex(this,side,index);
}

void TPZGeoElPi3d::NewMidSideNode(int side,int &index) {
	TPZRefPyramid::NewMidSideNode(this,side,index);
}

int TPZGeoElPi3d::SideDimension(int side) {
	return TPZShapePiram::SideDimension(side);
}

/*
TPZGeoElSide TPZGeoElPi3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a 4  targetdimension deve ser 1
//se side =  5 a 12  targetdimension deve ser 2
//se side = 18      targetdimension deve ser 3
  if( (side<0 || side>18) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElPi3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;
  //side = 0 a 18
  switch(targetdimension) {//=1,2
	  case 1:
       if(father->NSides() == 15) return TPZGeoElSide();//o pai é um tetraedro 
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,5);
       	 if(side==3) return TPZGeoElSide(this,8);
          if(side==4) return TPZGeoElSide(this,9);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,5);
       	 if(side==2) return TPZGeoElSide(this,6);
          if(side==4) return TPZGeoElSide(this,10);
       } else if(this == father->SubElement(2)) {
       	 if(side==1) return TPZGeoElSide(this,6);
       	 if(side==3) return TPZGeoElSide(this,7);
          if(side==4) return TPZGeoElSide(this,11);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(this,8);
       	 if(side==2) return TPZGeoElSide(this,7);
          if(side==4) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(4)) {
       	 if(side==0) return TPZGeoElSide(this,9);
       	 if(side==1) return TPZGeoElSide(this,10);
          if(side==2) return TPZGeoElSide(this,11);
          if(side==3) return TPZGeoElSide(this,12);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(13,14);
          if((side==1 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,17);
          if((side==3 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(14,17);
          if((side==4 || side==9) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==2 || side==6 || side==7) return TPZGeoElSide(this,13);
       	 if(side==10) return TPZGeoElSide(this,14);
       	 if(side==12) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(13,14);
          if((side==0 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,15);
          if((side==2 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(14,15);
          if((side==4 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          if(side==3 || side==7 || side==8) return TPZGeoElSide(this,13);
       	 if(side==9) return TPZGeoElSide(this,14);
       	 if(side==11) return TPZGeoElSide(this,15);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(13,15);
          if((side==1 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,16);
          if((side==3 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==4 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          if(side==0 || side==5 || side==8) return TPZGeoElSide(this,13);
       	 if(side==10) return TPZGeoElSide(this,15);
       	 if(side==12) return TPZGeoElSide(this,16);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(13,17);
          if((side==0 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,16);
          if((side==2 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==4 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          if(side==1 || side==5 || side==6) return TPZGeoElSide(this,13);
       	 if(side==9) return TPZGeoElSide(this,17);
       	 if(side==11) return TPZGeoElSide(this,16);
       } else if(this == father->SubElement(4)) {

          if(father->NSides() == 15) {//o pai é um tetraedro

          	 TPZGeoElT3d *sub = (TPZGeoElT3d *) father->SubElement(0);
             bestface = sub->BestDimensionSideOfTwoFaces(10,11);
             if(side==0 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(10,13);
             if(side==3 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(11,13);
             if(side==4 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==8)             return TPZGeoElSide(sub,10);
             sub = (TPZGeoElT3d *) father->SubElement(3);
             bestface = sub->BestDimensionSideOfTwoFaces(11,12);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(12,13);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==6)             return TPZGeoElSide(sub,12);
             if(side== 5 || side== 9 || side==10) return TPZGeoElSide(this,14);
             if(side== 7 || side==11 || side==12) return TPZGeoElSide(this,16);
             return TPZGeoElSide();
          }
          bestface = BestDimensionSideOfTwoFaces(14,17);
          if((side==0 || side==9) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(14,15);
          if((side==1 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==2 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==3 || side==12) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==5) return TPZGeoElSide(this,14);
       	 if(side==6) return TPZGeoElSide(this,15);
       	 if(side==7) return TPZGeoElSide(this,16);
       	 if(side==8) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(5)) {

          if(father->NSides() == 15) {//o pai é um tetraedro

          	 TPZGeoElT3d *sub = (TPZGeoElT3d *) father->SubElement(1);
             bestface = sub->BestDimensionSideOfTwoFaces(10,11);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(11,12);
             if(side==0 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(10,12);
             if(side==4 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==5)             return TPZGeoElSide(sub,11);
             sub = (TPZGeoElT3d *) father->SubElement(2);
             bestface = sub->BestDimensionSideOfTwoFaces(10,13);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(12,13);
             if(side==3 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==7)             return TPZGeoElSide(sub,13);
             if(side== 8 || side== 9 || side==12) return TPZGeoElSide(this,17);
             if(side== 6 || side==10 || side==11) return TPZGeoElSide(this,15);
             return TPZGeoElSide();
          }
          TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(0);
          bestface = sub->BestDimensionSideOfTwoFaces(14,17);
          if(side==1 && bestface==17) return TPZGeoElSide(father->SubElement(9),12);
          sub = (TPZGeoElPi3d *) father->SubElement(3);
          bestface = sub->BestDimensionSideOfTwoFaces(16,17);
          if(side==2 && bestface==17) return TPZGeoElSide(father->SubElement(9),12);
          sub = (TPZGeoElPi3d *) father->SubElement(2);
          bestface = sub->BestDimensionSideOfTwoFaces(15,16);
          if(side==3 && bestface==16) return TPZGeoElSide(father->SubElement(8),13);
       	 if(side==5) return TPZGeoElSide(father->SubElement(4),14);//(this,14)
       	 if(side==7) return TPZGeoElSide(father->SubElement(4),16);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,18);//0<=side<=18
  }//switch
  return TPZGeoElSide();
}
*/
void TPZGeoElPi3d::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides) {
	TPZStack<int> high;
	TPZShapePiram::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}
}

/*
int TPZGeoElPi3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

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
    if(levnum1 > levnum2) return face1;
    if(levnum2 > levnum1) return face2;
    return 0;
}
*/

void TPZGeoElPi3d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapePiram::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapePiram::SideConnectLocId(side,is));
}

void TPZGeoElPi3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
	TPZFMatrix xco(3,TPZGeoPyramid::NNodes);
	int i,j;
	for(i=0; i<TPZShapePiram::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoPyramid::Jacobian(xco,param,jacobian,axes,detjac,jacinv);
}

void TPZGeoElPi3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
	TPZFMatrix xco(3,TPZGeoPyramid::NNodes);
	int i,j;
	for(i=0; i<TPZShapePiram::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoPyramid::X(xco,loc,result);
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/


/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElPi3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
	TPZRefPyramid::Divide(this,SubElVec);
}


int TPZGeoElPi3d::NSubElements() {
	return TPZRefPyramid::NSubEl;
}

/*
TPZGeoElSide TPZGeoElPi3d::SideSubElement(int side,int position) {
   if (position<0 || position>10 || side <0 ||side>18) {
   	PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
   }                              //fSubEl[is]
   if(side==18) {
      if(position > 5) return TPZGeoElSide(SubElement(position),14);//centro : tetraedro
      return TPZGeoElSide(SubElement(position),18);//centro : pirâmides
   }
   if(side<5) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>4 && side<13) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
      	int s = side-5;
         return TPZGeoElSide(SubElement(TPZCompElPi3d::SideNodes[s][position]),side);
      }
   }
   if(side>12) {//faces
       if(position<0 || position>4) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
      	int s = side-13;//s = 0,1,2,3,4
         if(s==1 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),11);
         if(s==2 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),10);
         if(s==3 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),13);
         if(s==4 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),12);
         return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][position]),side);
      }
   }
   return TPZGeoElSide();
}

void TPZGeoElPi3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 18) {
      PZError << "TPZGeoElPi3d::SideSubElements called for side " << side << endl;
      return;
   }
   if(side==18) {
      sub.Resize(10);
      for(int i=0;i<10;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<5) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>4 && side<13) {//lados
      int s = side-5;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElPi3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElPi3d::SideNodes[s][1]];
      return;
   }
   if(side>12) {//faces
      int s = side-13;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElPi3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElPi3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElPi3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElPi3d::FaceSons[s][3]];
   }
}

TPZGeoElSide TPZGeoElPi3d::Father(int side) {

   if(!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i,nsides = fFather->NSides();
   if(nsides == 15)//o pai é um tetraedro
   	for(i=4;i<6;i++)  if(fFather->SubElement(i) == this) whichsub = i;
   if(nsides == 19)//o pai é uma pirâmide
   	for(i=0;i<6;i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {
	   PZError << "TPZGeoElPi3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   if(nsides == 15) {//face da pirâmide cujo pai é um tetraedro
      //if(side >13) {
         if(whichsub == 4 && (side==14 || side==16)) return TPZGeoElSide(fFather,side-3); else
         if(whichsub == 5 && (side==15 || side==17)) return TPZGeoElSide(fFather,side-5);
				 if(side==18) return TPZGeoElSide(fFather,14);//pai tetraedro para o interior do filho pirâmide
      //}
      return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 10
   //os filhos interiores não tém pai associados a seus cantos
   if((side<5 && side == whichsub) || side==18) return TPZGeoElSide(fFather,side);//cantos
   //lados
   if(whichsub == 0 && (side== 5 || side== 8 || side== 9)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side== 5 || side== 6 || side==10)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side== 6 || side== 7 || side==11)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side== 7 || side== 8 || side==12)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side== 9 || side==10 || side==11 || side==12)) return TPZGeoElSide(fFather,side);
   //faces
   if(whichsub == 0 && (side==13 || side==14 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==13 || side==14 || side==15)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==13 || side==15 || side==16)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side==13 || side==16 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side==14 || side==15 || side==16 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 &&  side==11)                          return TPZGeoElSide(fFather,11);
   if(whichsub == 7 &&  side==10)                          return TPZGeoElSide(fFather,10);
   if(whichsub == 8 &&  side==13)                          return TPZGeoElSide(fFather,13);
   if(whichsub == 9 &&  side==12)                          return TPZGeoElSide(fFather,12);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElPi3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 4
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElPi3d::GetSubElement subelement does not contain refnode" << endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
      int nnod = 4,son;
      for(son=6;son<10;son++) if(fSubEl[son] == subel) break;
      if(son == 10) nnod = 5;
		for(k = 0; k < refnode.NElements(); k++) {//se o subelemento k do thisside tiver o canto refnode[k]
		   for(j=0;j<nnod;j++) {  //este é um subelemento procurado
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   if(side > 13 && side < 18) sub[3] = SideSubElement(side,3);
   if(side == 18) {
   	sub[5] = SideSubElement(side,5);
      sub[6] = SideSubElement(side,6);
   	sub[7] = SideSubElement(side,7);
      sub[8] = SideSubElement(side,8);
      sub[9] = SideSubElement(side,9);
   }
   return;
}
*/
/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
//transforma piramide para piramide ou pirâmide para tetraedro
/*
void TPZGeoElPi3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform &t) {
   if(!fFather || side > 18) return;
   int whichsub = -1;
   int i,nsides = fFather->NSides();
   if(nsides == 15)//o pai é um tetraedro
   	for(i=4;i<6;i++)  if(fFather->SubElement(i) == this) whichsub = i;
   if(nsides == 19)//o pai é uma pirâmide
   	for(i=0;i<6;i++) if(fFather->SubElement(i) == this) whichsub = i;
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
   	cout << "TPZGeoElPi3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 18) {//pai mestre para filho
      mult(0,0) = 0.5;//ou transformação entre elemento mestre do filho
      mult(1,1) = 0.5;//para o elemento mestre do pai
      mult(2,2) = 0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
         case 0:
         	sum(0,0) = -.5;
            sum(1,0) = -.5;
            break;
         case 1:
         	sum(0,0) =  .5;
            sum(1,0) = -.5;
            break;
         case 2:
         	sum(0,0) = .5;
            sum(1,0) = .5;
            break;
         case 3:
         	sum(0,0) = -.5;
            sum(1,0) =  .5;
            break;
         case 4:
            if(nsides==15) {//pai tetraedro
               mult.Zero();
               mult(0,1) = -0.25;
               mult(0,2) = -0.25;
               mult(1,1) =  0.25;
               mult(1,2) = -0.25;
               mult(2,0) =  0.25;
               mult(2,2) =  0.25;
               sum(0,0)  = .25;
               sum(1,0)  = .25;
               sum(2,0)  = .25;
            } else {//pai pirâmide
               sum(2,0) = .5;
            }
         case 5:
         	if(nsides==19) {//pai pirâmide
               mult(0,0) *= -1.;
               mult(2,2) *= -1.;
               sum(2,0)   = .5;
            } else {//pai tetraedro
               mult.Zero();
               mult(0,1) = -0.25;
               mult(0,2) =  0.25;
               mult(1,1) =  0.25;
               mult(1,2) =  0.25;
               mult(2,0) = -0.25;
               mult(2,2) = -0.25;
               sum(0,0)  = .25;
               sum(1,0)  = .25;
               sum(2,0)  = .25;
            }
      }
   } else if(side>12) {//face do pai para face do filho
      whichsub = -1;
      int s = side-13;//lado da pirâmide
      if(nsides == 15) {//pai tetraedro
      	i = TPZCompElPi3d::MiddleFace[s-1]-10;//face do tetraedro que contém a face da pirâmide atual
         if(fFather->SubElement(TPZCompElT3d::FaceSons[i][3]) == this) whichsub = 3;
      }
      if(nsides == 19) {//pai pirâmide
         int n = 3;
         if(side==13) n = 4;
         for(i=0; i<n; i++) if(fFather->SubElement(TPZCompElPi3d::FaceSons[s][i]) == this) whichsub = i;
      }
      if(whichsub == -1) return;
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
   	if(side == 13) {//face quadrilateral
			sum(0,0) = -0.5;
         sum(1,0) = -0.5;
         switch(whichsub) {//o atual é a face numero whichsub do filho dentro da face do pai
            case 0:
               break;
            case 1:
               sum(0,0) = .5;
               break;
            case 2:
               sum(0,0) = .5;
               sum(1,0) = .5;
               break;
            case 3:
               sum(1,0) = .5;
         }
      } else {//face triangular
         switch(whichsub) {//o atual é o filho numero whichsub
            case 0://filho pirâmide e pai pirâmide casos 0,1,2
               break;
            case 1:
               sum(0,0) = 0.5;
               break;
            case 2:
               sum(1,0) = 0.5;
               break;
            case 3://filho piramide e pai tetraedro caso 3
               if(side==14 || side==16) {//basta com o side ja que o side do tetraedro acaba em 13 para as faces
                  mult(0,0) = 0.;
                  mult(0,1) =-0.5;
                  mult(1,0) = 0.5;
                  sum(0,0) = 0.5;
                  sum(1,0) = 0.;
               } else
               if(side==15) {
                  mult(0,0) =-0.5;
                  mult(1,0) = 0.5;
                  sum(0,0) = 0.5;
                  sum(1,0) = 0.;
               } else
               if(side==17) {
                  mult(0,1) = 0.5;
                  mult(1,1) =-0.5;
                  sum(0,0) = 0.;
                  sum(1,0) = 0.5;
               }
         }
      }
   } else if(side>4) {//e side < 13
      whichsub = -1;
      int s = side-5;//0 a 7
      for(i=0; i<2; i++) {
         if(nsides==15) if(fFather->SubElement(TPZCompElT3d::SideNodes[s-1][i]) == this) whichsub = i;
         if(nsides==19) if(fFather->SubElement(TPZCompElPi3d::SideNodes[s][i]) == this) whichsub = i;
      }
      if(whichsub == -1) return;
   	mult(0,0) = 0.5;
      if(whichsub==0) sum(0,0) = -0.5;//subelemento 0 do lado
      if(whichsub==1) sum(0,0) =  0.5;//subelemento 1 do lado
   }
   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform(side,father,t);
}
*/
TPZTransform TPZGeoElPi3d::SideToSideTransform(int sidefrom,int sideto) {
	return TPZShapePiram::SideToSideTransform(sidefrom,sideto);
}

TPZGeoEl *TPZGeoElPi3d::CreateBCGeoEl(int side,int bc) {
	TPZGeoEl *gel = TPZGeoPyramid::CreateBCGeoEl(this,side,bc);
	return gel;
}

/*
void TPZGeoElPi3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(4,-1);
   if((face>-1 && face<5) || (face>12 && face<18)) {
   	if(face>12) face = face-13;
      ids[0] = NodeIndex(TPZCompElPi3d::FaceNodes[face][0]);
      ids[1] = NodeIndex(TPZCompElPi3d::FaceNodes[face][1]);
      ids[2] = NodeIndex(TPZCompElPi3d::FaceNodes[face][2]);
		ids[3] = NodeIndex(TPZCompElPi3d::FaceNodes[face][3]);
      if(face > 0) ids.Resize(3);
      return;
   }
 	cout << "TPZCompElPi3d::NodeFaceIds bad side , side = " << face << endl;
}
*/

//cada lado do filho {s0,s1,s2,...,s19} está contido num lado do pai: sn=fatside

TPZGeoElSide TPZGeoElPi3d::Father2(int side){
	if (!fFather ) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0 || !fFather) return TPZGeoElSide();
	int fathsid = fFather->FatherSide(side,son);
	return TPZGeoElSide(fFather,fathsid);
}

int TPZGeoElPi3d::FatherSide(int side, int son){
	return TPZRefPyramid::FatherSide(side,son);
}

void TPZGeoElPi3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
	TPZRefPyramid::GetSubElements(this,side,subel);
}

int TPZGeoElPi3d::NSideSubElements2(int side) {
	return TPZRefPyramid::NSideSubElements(side);
}


TPZTransform TPZGeoElPi3d::BuildTransform2(int side, TPZGeoEl *father, TPZTransform &t){//Augusto:09/01/01

	if(side<0 || side>18 || !father){
  	PZError << "TPZGeoElPi3d::BuildTransform2 side out of range or father null\n";
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


REAL TPZGeoElPi3d::MidSideNode[19][3] = {
/*00*/{-1.,-1.},   /*01*/{1.,-1.},   /*02*/{1.,1.},/*03*/{-1.,1.},/*04*/{0.,0.,1.},
/*05*/{ 0.,-1.},   /*06*/{1., 0.},   /*07*/{0.,1.},/*08*/{-1.,0.},
/*09*/{-.5,-.5,.5},/*10*/{.5,-.5,.5},/*11*/{.5,.5,.5},/*12*/{-.5,.5,.5},
/*13*/{0.,  0. ,  0. },/*14*/{  0.  ,-2./3.,1./3.},/*15*/{2./3.,0.,1./3.},
/*16*/{0.,2./3.,1./3.},/*17*/{-2./3.,  0.  ,1./3.},/*18*/{  0. ,0.,1./5.} };


int TPZGeoElPi3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<10;sn++){
    TPZGeoEl *son = subs[sn];
    int nsides = son->NSides();
    for(sd=0;sd<nsides;sd++){
      TPZTransform telsd(0,0);
      if(nsides==15){
        ps[0] = TPZGeoElT3d::MidSideNode[sd][0];//tetraedro
        ps[1] = TPZGeoElT3d::MidSideNode[sd][1];
        ps[2] = TPZGeoElT3d::MidSideNode[sd][2];
        telsd = TPZShapeTetra::TransformElementToSide(sd);
      } else if(nsides==19){
        ps[0] = MidSideNode[sd][0];//pirâmide
        ps[1] = MidSideNode[sd][1];
        ps[2] = MidSideNode[sd][2];
        telsd = TPZShapePiram::TransformElementToSide(sd);
      }
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform tson(son->SideDimension(sd));
      TPZTransform t = gel->BuildTransform2(sd,son,tson);//para não furar pirâmide com pai tetraedro
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = TPZShapePiram::TransformSideToElement(sdfat); //else
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(nsides-1).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao errada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(nsides-1).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(nsides-1).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}

void TPZGeoElPi3d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >9){
    PZError << "TPZGeoElPi3d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}
   


TPZIntPoints * TPZGeoElPi3d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoPyramid::CreateSideIntegrationRule(side,order);
}

TPZTransform TPZGeoElPi3d::GetTransform(int side,int son) {
	return TPZRefPyramid::GetTransform(side,son);
}

void TPZGeoElPi3d::CenterPoint(int side, TPZVec<REAL> &masscent){

  if(side < 0 || side > (NSides()-1)){
    PZError << "TPZGeoElPi3d::CenterPoint error side = " << side << endl;
    return;
  }
  TPZShapePiram::CenterPoint(side,masscent);
}

REAL TPZGeoElPi3d::RefElVolume(){
  return TPZShapePiram::RefElVolume();
}
