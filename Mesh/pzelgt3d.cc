// -*- c++ -*-
//METHODS DEFINITION FOR CLASS ELEMGTD
#include "pzelgt3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgt2d.h"
#include "pzelgpi3d.h"
#include "pzelct3d.h"
#include "pzelcpi3d.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElT3d *gel,TPZCompMesh &mesh,int &index) {
  return new TPZCompElT3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElT3d::fp)(TPZGeoElT3d *,TPZCompMesh &,int &) = CreateEl;

TPZGeoElT3d::TPZGeoElT3d(int id,TPZVec<int> &nodeindexes,int matid,TPZGeoMesh &mesh):
  TPZGeoEl(id,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElT3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<6;i++) fSubEl[i] = 0;
}

TPZGeoElT3d::TPZGeoElT3d(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElT3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<6;i++) fSubEl[i] = 0;
}

TPZGeoElT3d::~TPZGeoElT3d() {}

void TPZGeoElT3d::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
	TPZGeoTetrahedra::Shape(pt,phi,dphi);
}

TPZGeoElT3d *TPZGeoElT3d::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElT3d(nodeindexes,matind,mesh);
}

int TPZGeoElT3d::NNodes() {
	return TPZShapeTetra::NNodes;
}
int TPZGeoElT3d::NodeIndex(int node) {
  if(node<0 || node>3) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElT3d::NSideNodes(int side) {
	return TPZShapeTetra::NSideNodes(side);
}

int TPZGeoElT3d::SideNodeIndex(int side,int node) {
	int loc = TPZShapeTetra::SideNodeLocId(side,node);
	if(loc<0) return loc;
	return fNodeIndexes[loc];
}

void TPZGeoElT3d::MidSideNodeIndex(int side,int &index) {
	TPZRefTetrahedra::MidSideNodeIndex(this,side,index);
}

void TPZGeoElT3d::NewMidSideNode(int side,int &index) {
	TPZRefTetrahedra::NewMidSideNode(this,side,index);
}

int TPZGeoElT3d::SideDimension(int side) {
	return TPZShapeTetra::SideDimension(side);
}
/*
TPZGeoElSide TPZGeoElT3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a 3  targetdimension deve ser 1
//se side =  4 a 9  targetdimension deve ser 2
//se side = 14      targetdimension deve ser 3
  if( (side<0 || side>14) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElT3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;
  //side = 0 a 14
  switch(targetdimension) {//=1,2
	  case 1:
       if(father->NSides() == 19) {//o pai é uma pirâmide
          if(this == father->SubElement(6)) {
          	 if(side==0) return TPZGeoElSide(father->SubElement(4), 9);
             if(side==1) return TPZGeoElSide(father->SubElement(0), 5);
             if(side==2) return TPZGeoElSide(this,5);//canto do mesmo id
          } else if(this == father->SubElement(7)) {
             if(side==1) return TPZGeoElSide(father->SubElement(4),10);
             if(side==2) return TPZGeoElSide(father->SubElement(4),11);
             if(side==3) return TPZGeoElSide(this,7);

          } else if(this == father->SubElement(8)) {
             if(side==1) return TPZGeoElSide(this,5);
             if(side==2) return TPZGeoElSide(father->SubElement(3), 7);

          } else if(this == father->SubElement(9)) {
          	 if(side==0) return TPZGeoElSide(this,7);
             if(side==3) return TPZGeoElSide(father->SubElement(0),8);
          }
       	return TPZGeoElSide();
       }
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,4);
       	 if(side==2) return TPZGeoElSide(this,6);
          if(side==3) return TPZGeoElSide(this,7);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,4);
       	 if(side==2) return TPZGeoElSide(this,5);
          if(side==3) return TPZGeoElSide(this,8);
       } else if(this == father->SubElement(2)) {
       	 if(side==0) return TPZGeoElSide(this,6);
       	 if(side==1) return TPZGeoElSide(this,5);
          if(side==3) return TPZGeoElSide(this,9);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(this,7);
       	 if(side==1) return TPZGeoElSide(this,8);
          if(side==2) return TPZGeoElSide(this,9);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(10,11);
          if((side==1 || side==4) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(10,13);
          if((side==2 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(11,13);
          if((side==3 || side==7) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 5) return TPZGeoElSide(this,10);
       	 if(side== 8) return TPZGeoElSide(this,11);
       	 if(side== 9) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(10,11);
          if((side==0 || side==4) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(10,12);
          if((side==2 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(11,12);
          if((side==3 || side==8) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 6) return TPZGeoElSide(this,10);
       	 if(side== 7) return TPZGeoElSide(this,11);
       	 if(side== 9) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(10,12);
          if((side==1 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(10,13);
          if((side==0 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(12,13);
          if((side==2 || side==9) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 4) return TPZGeoElSide(this,10);
       	 if(side== 7) return TPZGeoElSide(this,13);
       	 if(side== 8) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(11,13);
          if((side==0 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(11,12);
          if((side==1 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(12,13);
          if((side==2 || side==9) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 4) return TPZGeoElSide(this,11);
       	 if(side== 5) return TPZGeoElSide(this,12);
       	 if(side== 6) return TPZGeoElSide(this,13);
       }
       //o pai é uma pirâmide
       if(father->NSides() == 19) {
          if(this == father->SubElement(6)) {

             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(4);//pois primeior subira pelo canto para aresta e logo de aresta para face
             bestface = sub->BestDimensionSideOfTwoFaces(14,17);
             if(side==0 && bestface) return TPZGeoElSide(sub,bestface);
             sub = (TPZGeoElPi3d *) father->SubElement(0);
             bestface = sub->BestDimensionSideOfTwoFaces(13,14);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);//para side=5 a aresta do sub 0 nao é 5, entao nao pode, outro irmao o considerara
             if(side==2) return TPZGeoElSide(sub,13);//para side 5 as faces do atual que contem este side
             if(side==4 || side==7 || side==8) return TPZGeoElSide(this,11);

          } else if(this == father->SubElement(7)) {
             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(4);//1
             bestface = sub->BestDimensionSideOfTwoFaces(14,15);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(15,16);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             sub = (TPZGeoElPi3d *) father->SubElement(1);
             if(side==3 || side==7) return TPZGeoElSide(sub,13);
             if(side==4 || side==5 || side==6) return TPZGeoElSide(this,10);

          } else if(this == father->SubElement(8)) {

             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(3);//2
             if(side==1) return TPZGeoElSide(sub,13);
             bestface = sub->BestDimensionSideOfTwoFaces(13,16);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==6 || side==7 || side==9) return TPZGeoElSide(this,13);
          } else if(this == father->SubElement(9)) {

             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(0);
             bestface = sub->BestDimensionSideOfTwoFaces(13,17);
             if(side==3 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==7) return TPZGeoElSide(sub,13);
             if(side==5 || side==8 || side==9) return TPZGeoElSide(this,12);
          }
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,14);//0<=side<=14
  }//switch
  return TPZGeoElSide();
}
*/
void TPZGeoElT3d::AllHigherDimensionSides(int side,int targetdimension, TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> high;
	TPZShapeTetra::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}
}

/*
int TPZGeoElT3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

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

void TPZGeoElT3d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapeTetra::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapeTetra::SideConnectLocId(side,is));

}

void TPZGeoElT3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
	TPZFMatrix xco(3,TPZGeoTetrahedra::NNodes);
	int i,j;
	for(i=0; i<TPZShapeTetra::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoTetrahedra::Jacobian(xco,param,jacobian,axes,detjac,jacinv);
}

void TPZGeoElT3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
	TPZFMatrix xco(3,TPZGeoTetrahedra::NNodes);
	int i,j;
	for(i=0; i<TPZShapeTetra::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoTetrahedra::X(xco,loc,result);
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/


/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElT3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
	TPZRefTetrahedra::Divide(this,SubElVec);
}


int TPZGeoElT3d::NSubElements() {
	return TPZRefTetrahedra::NSubEl;
}

/*
TPZGeoElSide TPZGeoElT3d::SideSubElement(int side,int position) {
   if (position<0 || position>6 || side <0 ||side>14) {
   	PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
   }                              //fSubEl[is]
   if(side==14) {
      if(position > 3) return TPZGeoElSide(SubElement(position),18);//centro : pirâmides
      return TPZGeoElSide(SubElement(position),14);//centro : tetraedros
	}
   if(side<4) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>3 && side<10) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
      	int s = side-4;
         return TPZGeoElSide(SubElement(TPZCompElT3d::SideNodes[s][position]),side);
      }
   }
   if(side>9) {//faces
       if(position<0 || position>3) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
      	int s = side-10;//s = 0,1,2,3
         if(s==0 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),15);//lado do subelemento
         if(s==1 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),14);//interior à face
         if(s==2 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),17);//nestes casos uma
         if(s==3 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),16);//piramide
         return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][position]),side);
      }
   }
   return TPZGeoElSide();
}

void TPZGeoElT3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 14) {
      PZError << "TPZGeoElT3d::SideSubElements called for side " << side << endl;
      return;
   }
   if(side==14) {
      sub.Resize(6);
      for(int i=0;i<6;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<4) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>3 && side<10) {//lados
      int s = side-4;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElT3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElT3d::SideNodes[s][1]];
      return;
   }
   if(side>9) {//faces
      int s = side-10;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElT3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElT3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElT3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElT3d::FaceSons[s][3]];
   }
}

TPZGeoElSide TPZGeoElT3d::Father(int side) {

   if(!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i,nsubel = 6;
   if(fFather->NSides() == 19) nsubel = 10;//pai é pirâmide
   for(i=0; i<nsubel; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {
	   PZError << "TPZGeoElT3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   if(fFather->NSides() == 19) {//o pai é uma pirâmide
      if(whichsub == 6 &&  side==11) return TPZGeoElSide(fFather,14);
      if(whichsub == 7 &&  side==10) return TPZGeoElSide(fFather,15);
      if(whichsub == 8 &&  side==13) return TPZGeoElSide(fFather,16);
      if(whichsub == 9 &&  side==12) return TPZGeoElSide(fFather,17);
      if(side==14) return TPZGeoElSide(fFather,18);//pai do tetraedro pelo interior
      return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 6
   //os filhos interiores não tém pai associados a seus cantos
   if((side<4 && side == whichsub) || side==14) return TPZGeoElSide(fFather,side);//cantos
   //lados
   if(whichsub == 0 && (side== 4 || side== 6 || side== 7)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side== 4 || side== 5 || side== 8)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side== 5 || side== 6 || side== 9)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side== 7 || side== 8 || side== 9)) return TPZGeoElSide(fFather,side);
   //faces
   if(whichsub == 0 && (side==10 || side==11 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==10 || side==11 || side==12)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==10 || side==12 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side==11 || side==12 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side==14 || side==16))             return TPZGeoElSide(fFather,side-3);
   if(whichsub == 5 && (side==15 || side==17))             return TPZGeoElSide(fFather,side-5);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElT3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 3
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElT3d::GetSubElement subelement does not contain refnode" << endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k = 0; k < refnode.NElements(); k++) {//k<nsub?
		   for(j=0;j<4;j++) {//4 vantos do pai
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);//sub[k]?
            }
         }
      }
   }
   if(side > 9 && side < 14) sub[3] = SideSubElement(side,3);
   if(side == 14) {
   	sub[4] = SideSubElement(side,4);
      sub[5] = SideSubElement(side,5);
   }
   return;
}
*/
/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
//transforma tetraedro para tetraedro ou tetraedro para pirâmide
//subelementos 0,1,2,3 do tetraedro -> tetraedro
//subelementos 6,7,8,9 da tetraedro -> pirâmide
//transformação entre lado filho e lado pai, o side é o side do filho
//o lado do filho esta contido no lado do pai
/*
void TPZGeoElT3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform &t) {
   if(!fFather || side > 14) return;
   int whichsub = -1;
   int i,nsides = fFather->NSides();
   if(nsides == 15)//o pai é um tetraedro
   	for(i=0;i<4;i++)  if(fFather->SubElement(i) == this) whichsub = i;
   if(nsides == 19)//o pai é uma pirâmide
   	for(i=6;i<10;i++) if(fFather->SubElement(i) == this) whichsub = i;
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
   	cout << "TPZGeoElT3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side==14) {//pai para filho ou filho mestre para pai mestre
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      mult(2,2) = 0.5;
      switch(whichsub) {//o atual é o tetraedro filho numero whichsub
         case 0:
            break;
         case 1:
            sum(0,0) = .5;
            break;
         case 2:
            sum(1,0) = .5;
            break;
         case 3:
         	sum(2,0) = .5;
            break;
         case 6:
         	mult.Zero();
            mult(0,0) =  0.5;
            mult(0,1) =  1.;
            mult(0,2) =  0.5;
            mult(1,2) =  0.5;
            mult(2,0) = -0.5;
            mult(2,2) = -0.5;
         	sum(0,0)  = -.5;
            sum(1,0)  = -.5;
            sum(2,0)  =  .5;
            break;
         case 7:
         	mult.Zero();
            mult(0,0) = -0.5;
            mult(0,1) = -0.5;
            mult(0,2) = -1.;
            mult(1,0) = -0.5;
            mult(1,1) =  0.5;
            mult(2,0) =  0.5;
            mult(2,1) =  0.5;
         	sum(0,0)  = 1.;
            break;
         case 8:
         	mult.Zero();
            mult(0,0) =  0.5;
            mult(0,1) =  0.5;
            mult(0,2) =  1.;
            mult(1,0) = -0.5;
            mult(1,1) =  0.5;
            mult(2,0) = -0.5;
            mult(2,1) = -0.5;
         	sum(0,0)  = -.5;
            sum(1,0)  =  .5;
            sum(2,0)  =  .5;
            break;
         case 9:
         	mult.Zero();
            mult(0,0) = -0.5;
            mult(0,1) = -0.5;
            mult(0,2) = -1.;
            mult(1,0) = -0.5;
            mult(1,1) =  0.5;
            mult(2,0) =  0.5;
            mult(2,1) =  0.5;
      }
   } else if(side>9) {//face do pai para face do filho
         whichsub = -1;
         int s = side-10;;
         if(nsides == 15) {
            for(i=0; i<3; i++) if(fFather->SubElement(TPZCompElT3d::FaceSons[s][i]) == this) whichsub = i;
         }
         if(nsides == 19) {
               i = TPZCompElT3d::MiddleFace[s]-13;//face da pirâmide que contem a face do tetraedro atual
            	if(fFather->SubElement(TPZCompElPi3d::FaceSons[i][3]) == this) whichsub = 3;
         }
         if(whichsub == -1) return;
         mult(0,0) = 0.5;
         mult(1,1) = 0.5;
         switch(whichsub) {//o atual é o filho numero whichsub
            case 0:        //ambos são tetraedros para os casos 0,1,2
               break;
            case 1:
               sum(0,0) = 0.5;
               break;
            case 2:
               sum(1,0) = 0.5;
               break;
            case 3://filho tetraedro e pai pirâmide
               if(side==11 || side==13) {//basta com o side do filho já
                  mult(0,1) = 0.5;        //que o pai é pirâmide
                  mult(1,0) =-0.5;
                  mult(1,1) = 0.;
                  sum(1,0)  = 0.5;
               } else
               if(side==10) {
                  mult(0,0) =-0.5;
                  mult(1,0) = 0.5;
                  sum(0,0)  = 0.5;
               } else
               if(side==12) {
                  mult(0,1) = 0.5;
                  mult(1,1) =-0.5;
                  sum(1,0)  = 0.5;
               }
         }
   } else if(side>3) {//4 a 8
      whichsub = -1;
      int s = side-4;
      for(i=0; i<2; i++) {
         if(nsides==15) if(fFather->SubElement(TPZCompElT3d::SideNodes[s][i]) == this) whichsub = i;
         if(nsides==19) if(fFather->SubElement(TPZCompElPi3d::SideNodes[s+1][i]) == this) whichsub = i;
      }
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
TPZTransform TPZGeoElT3d::SideToSideTransform(int sidefrom,int sideto) {
	return TPZShapeTetra::SideToSideTransform(sidefrom,sideto);
}

TPZGeoEl *TPZGeoElT3d::CreateBCGeoEl(int side,int bc) {
	TPZGeoEl *gel = TPZGeoTetrahedra::CreateBCGeoEl(this,side,bc);
	return gel;
}

/*
void TPZGeoElT3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(3,-1);
   if((face>-1 && face<4) || (face>9 && face<14)) {
   	if(face>9) face = face-10;
      ids[0] = NodeIndex(TPZCompElT3d::FaceNodes[face][0]);
      ids[1] = NodeIndex(TPZCompElT3d::FaceNodes[face][1]);
      ids[2] = NodeIndex(TPZCompElT3d::FaceNodes[face][2]);
      return;
   }
 	cout << "TPZCompElT3d::NodeFaceIds bad side , side = " << face << endl;
}
*/
//cada lado do filho esta contido em que lado do pai?
static int fatherside[6][19] = {
/*00*/{0,4,6,7,4,10,6,7,11,13,10,11,14,13,14,-1,-1,-1,-1},
/*01*/{4,1,5,8,4,5,10,11,8,12,10,11,12,14,14,-1,-1,-1,-1},
/*02*/{6,5,2,9,10,5,6,13,12,9,10,14,12,13,14,-1,-1,-1,-1},
/*03*/{7,8,9,3,11,12,13,7,8,9,14,11,12,13,14,-1,-1,-1,-1},
/*04*/{4,8,9,6,7,11,12,13,10,11,11,13,13,14,11,14,13,14,14},
/*05*/{8,4,6,9,5,11,10,13,12,12,10,10,12,14,14,10,14,12,14},
};
static int fatherside2[4][19] = {//tetraedro com pai pirâmide
/*06*/{9,5,13,10,14,13,18,14,14,18,18,14,18,18,18,-1,-1,-1,-1},
/*07*/{6,10,11,13,15,15,15,13,18,18,15,18,18,18,18,-1,-1,-1,-1},
/*08*/{12,13,7,11,18,13,16,16,18,16,18,18,18,16,18,-1,-1,-1,-1},
/*09*/{13,9,12,8,18,17,18,13,17,17,18,18,17,18,18,-1,-1,-1,-1} };

TPZGeoElSide TPZGeoElT3d::Father2(int side){
	if (!fFather ) return TPZGeoElSide();	
	int son = WhichSubel();
	if(son<0) return TPZGeoElSide();
	int fathsid = fFather->FatherSide(side,son);
	return TPZGeoElSide(fFather,fathsid);
}

int TPZGeoElT3d::FatherSide(int side, int son){
	return TPZRefTetrahedra::FatherSide(side,son);
}

void TPZGeoElT3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

	TPZRefTetrahedra::GetSubElements(this,side,subel);
}

int TPZGeoElT3d::NSideSubElements2(int side) {
	return TPZRefTetrahedra::NSideSubElements(side);
}


TPZTransform TPZGeoElT3d::BuildTransform2(int side, TPZGeoEl *father, TPZTransform &t){//Augusto:09/01/01

	if(side<0 || side>18 || !father){
  	PZError << "TPZGeoElT3d::BuildTransform2 side out of range or father null\n";
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

REAL TPZGeoElT3d::MidSideNode[15][3] = {
/*00*/{.0,.0},/*01*/{1.,.0},/*02*/{0.,1.,.0},/*03*/{.0,0.,1.0},/*04*/{.5,.0,.0},
/*05*/{.5,.5},/*06*/{0.,.5},/*07*/{0.,0.,.5},/*08*/{.5,0.,0.5},/*09*/{.0,.5,.5},
/*10*/{1./3.,1./3., 0.  }  ,/*11*/{1./3., .0  ,1./3.},
/*12*/{1./3.,1./3.,1./3.}  ,/*13*/{ 0.  ,1./3.,1./3.},/*14*/{1./4.,1./4.,1./4.} };

int TPZGeoElT3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<6;sn++){
    TPZGeoEl *son = subs[sn];
    int nsides = son->NSides();
    for(sd=0;sd<nsides;sd++){
      TPZTransform telsd(0,0);
      if(nsides==15){
        ps[0] = MidSideNode[sd][0];//tetraedro
        ps[1] = MidSideNode[sd][1];
        ps[2] = MidSideNode[sd][2];
        telsd = TPZShapeTetra::TransformElementToSide(sd);
      } else if(nsides==19){
        ps[0] = TPZGeoElPi3d::MidSideNode[sd][0];//pirâmide
        ps[1] = TPZGeoElPi3d::MidSideNode[sd][1];
        ps[2] = TPZGeoElPi3d::MidSideNode[sd][2];
        telsd = TPZShapePiram::TransformElementToSide(sd);
      }
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform sont(son->SideDimension(sd));
      TPZTransform t = gel->BuildTransform2(sd,son,t);//para não furar pirâmide com pai tetraedro
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      /*if(nsides==15)*/
      telsd = TPZShapeTetra::TransformSideToElement(sdfat); //else
      //if(nsides==19) telsd = TPZShapePiram::TransformSideToElement(sdfat);
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

void TPZGeoElT3d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >5){
    PZError << "TPZGeoElT3d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}
 

TPZIntPoints * TPZGeoElT3d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoTetrahedra::CreateSideIntegrationRule(side,order);
}

TPZTransform TPZGeoElT3d::GetTransform(int side,int son) {
	return TPZRefTetrahedra::GetTransform(side,son);
}

void TPZGeoElT3d::CenterPoint(int side, TPZVec<REAL> &masscent){

  TPZShapeTetra::CenterPoint(side,masscent);
}

REAL TPZGeoElT3d::RefElVolume(){
  return TPZShapeTetra::RefElVolume();
}
