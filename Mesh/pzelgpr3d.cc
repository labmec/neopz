// -*- c++ -*-
//METHODS DEFINITION FOR CLASS ELEMPR3D
#include "pzelgpr3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelcpr3d.h"
#include "pzshapeprism.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElPr3d *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElPr3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElPr3d::fp)(TPZGeoElPr3d *,TPZCompMesh &,int &) = CreateEl;

TPZGeoElPr3d::TPZGeoElPr3d(int id,TPZVec<int> &nodeindexes,int matid,TPZGeoMesh &mesh):
  TPZGeoEl(id,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=6) {
    PZError << "TPZGeoElPr3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<6;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElPr3d::TPZGeoElPr3d(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=6) {
    PZError << "TPZGeoElPr3d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<6;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElPr3d::~TPZGeoElPr3d() {}

void TPZGeoElPr3d::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {

	TPZGeoPrism::Shape(pt,phi,dphi);
}

TPZGeoElPr3d *TPZGeoElPr3d::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElPr3d(nodeindexes,matind,mesh);
}

int TPZGeoElPr3d::NNodes() {
	return TPZShapePrism::NNodes;
}
int TPZGeoElPr3d::NodeIndex(int node) {
  if(node<0 || node>6) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElPr3d::NSideNodes(int side) {
	return TPZShapePrism::NSideNodes(side);
}

int TPZGeoElPr3d::SideNodeIndex(int side,int node) {
	int loc = TPZShapePrism::SideNodeLocId(side,node);
	if(loc<0) return loc;
	return fNodeIndexes[loc];
}

int TPZGeoElPr3d::SideNodeLocIndex(int side,int node) {
	return TPZShapePrism::SideNodeLocId(side,node);
}

void TPZGeoElPr3d::MidSideNodeIndex(int side,int &index) {
	TPZRefPrism::MidSideNodeIndex(this,side,index);
}

void TPZGeoElPr3d::NewMidSideNode(int side,int &index) {
	TPZRefPrism::NewMidSideNode(this,side,index);
}

int TPZGeoElPr3d::SideDimension(int side) {
	return TPZShapePrism::SideDimension(side);
}
/*
TPZGeoElSide TPZGeoElPr3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a 4  targetdimension deve ser 1
//se side =  5 a 12  targetdimension deve ser 2
//se side = 18      targetdimension deve ser 3
  if( (side<0 || side>20) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElPr3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;
  //side = 0 a 18
  switch(targetdimension) {//=1,2
	  case 1:
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,6);
       	 if(side==2) return TPZGeoElSide(this,8);
          if(side==3) return TPZGeoElSide(this,9);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,6);
       	 if(side==2) return TPZGeoElSide(this,7);
          if(side==4) return TPZGeoElSide(this,10);
       } else if(this == father->SubElement(2)) {
       	 if(side==0) return TPZGeoElSide(this,8);
       	 if(side==1) return TPZGeoElSide(this,7);
          if(side==5) return TPZGeoElSide(this,11);
       } else if(this == father->SubElement(3)) {
       	 if(side==3) return TPZGeoElSide(father->SubElement(0),8);
       	 if(side==5) return TPZGeoElSide(father->SubElement(0),6);
       	 if(side==4) return TPZGeoElSide(father->SubElement(1),7);
       } else if(this == father->SubElement(4)) {
       	 if(side==0) return TPZGeoElSide(this,9);
       	 if(side==4) return TPZGeoElSide(this,12);
          if(side==5) return TPZGeoElSide(this,14);
       } else if(this == father->SubElement(5)) {
       	 if(side==1) return TPZGeoElSide(this,10);
          if(side==3) return TPZGeoElSide(this,12);
       	 if(side==5) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(6)) {
       	 if(side==2) return TPZGeoElSide(this,11);
       	 if(side==3) return TPZGeoElSide(this,14);
          if(side==4) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(7)) {
       	 if(side==0) return TPZGeoElSide(father->SubElement(4),14);
       	 if(side==1) return TPZGeoElSide(father->SubElement(5),13);
       	 if(side==2) return TPZGeoElSide(father->SubElement(4),12);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==1 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,18);
          if((side==2 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,18);
          if((side==3 || side==9) && bestface) return TPZGeoElSide(this,bestface);
          if(side==7) return TPZGeoElSide(this,15);
       	 if(side==4 || side==10 || side==12) return TPZGeoElSide(this,16);
          if(side==5 || side==11 || side==14) return TPZGeoElSide(this,18);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==0 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,17);
          if((side==2 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==4 || side==10) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==8) return TPZGeoElSide(this,15);
          if(side==3 || side== 9 || side==12) return TPZGeoElSide(this,16);
          if(side==5 || side==11 || side==13) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(15,18);
          if((side==0 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,17);
          if((side==1 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(17,18);
          if((side==5 || side==11) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==6) return TPZGeoElSide(this,15);
          if(side==3 || side== 9 || side==14) return TPZGeoElSide(this,18);
          if(side==4 || side==10 || side==13) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(15,18);
          if(side==3 && bestface) return TPZGeoElSide(father->SubElement(0),bestface);
          bestface = BestDimensionSideOfTwoFaces(15,17);
          if(side==4 && bestface) return TPZGeoElSide(father->SubElement(1),bestface);
          bestface = BestDimensionSideOfTwoFaces(15,16);
          if(side==5 && bestface) return TPZGeoElSide(father->SubElement(0),bestface);
          if(side==12 || side==13 || side==14) return TPZGeoElSide(this,19);
          if(side== 0 || side== 9) return TPZGeoElSide(father->SubElement(6),18);//sub 6 dado que em ambos o canto é o zero
          if(side== 1 || side==10) return TPZGeoElSide(father->SubElement(5),17);//o canto é o mesmo
          //if(side== 2 || side==11) return TPZGeoElSide(father->SubElement(1),16);
          //como nao existe um irmao com o mesmo canto => nao existe HigherDimensioSide para este caso
       } else if(this == father->SubElement(4)) {

          bestface = BestDimensionSideOfTwoFaces(16,18);
          if((side==0 || side==9) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,19);
          if((side==4 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(18,19);
          if((side==5 || side==14) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==13) return TPZGeoElSide(this,19);
          if(side==1 || side== 6 || side==10) return TPZGeoElSide(this,16);
          if(side==2 || side== 8 || side==11) return TPZGeoElSide(this,18);
       } else if(this == father->SubElement(5)) {

          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==1 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,19);
          if((side==3 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(17,19);
          if((side==5 || side==13) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==14) return TPZGeoElSide(this,19);
          if(side==0 || side== 6 || side== 9) return TPZGeoElSide(this,16);
          if(side==2 || side== 7 || side==11) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(6)) {

          bestface = BestDimensionSideOfTwoFaces(17,18);
          if((side==2 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(18,19);
          if((side==3 || side==14) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(17,19);
          if((side==4 || side==13) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==12) return TPZGeoElSide(this,19);
          if(side==0 || side== 8 || side== 9) return TPZGeoElSide(this,18);
          if(side==1 || side== 7 || side==10) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(7)) {

          bestface = BestDimensionSideOfTwoFaces(18,19);
          if(side==0 && bestface) return TPZGeoElSide(father->SubElement(4),bestface);
          bestface = BestDimensionSideOfTwoFaces(17,19);
          if(side==1 && bestface) return TPZGeoElSide(father->SubElement(5),bestface);
          bestface = BestDimensionSideOfTwoFaces(16,19);
          if(side==2 && bestface) return TPZGeoElSide(father->SubElement(4),bestface);
          if(side==6 || side==7 || side==8) return TPZGeoElSide(this,15);
          if(side== 3 || side== 9) return TPZGeoElSide(father->SubElement(2),18);//o canto é o mesmo
          //if(side== 5 || side==11) return TPZGeoElSide(father->SubElement(4),16);//para este nao existe um irmao com o mesmo canto logo a transformacao deve ser consertada para este caso
          if(side== 4 || side==10) return TPZGeoElSide(father->SubElement(1),17);//o canto é o mesmo
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,20);//0<=side<=20
  }//switch
  return TPZGeoElSide();
}
*/
void TPZGeoElPr3d::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides) {
	TPZStack<int> high;
	TPZShapePrism::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}
}

/*
int TPZGeoElPr3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

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
//feito
void TPZGeoElPr3d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapePrism::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapePrism::SideConnectLocId(side,is));
}

//void TPZGeoElPr3d::SideMasterCo(int /*side*/,TPZFMatrix &/*coord*/) {
//}

//feito
void TPZGeoElPr3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
	TPZFMatrix xco(3,TPZGeoPrism::NNodes);
	int i,j;
	for(i=0; i<TPZShapePrism::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoPrism::Jacobian(xco,param,jacobian,axes,detjac,jacinv);
}

void TPZGeoElPr3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
	TPZFMatrix xco(3,TPZGeoPrism::NNodes);
	int i,j;
	for(i=0; i<TPZShapePrism::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoPrism::X(xco,loc,result);
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/


/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElPr3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
	TPZRefPrism::Divide(this,SubElVec);
}


int TPZGeoElPr3d::NSubElements() {
	return TPZRefPrism::NSubEl;
}

//int TPZGeoElPr3d::NSideSubElements(int side) {
//	return TPZRefPrism::NSideSubElements(side);
//}

/*
TPZGeoElSide TPZGeoElPr3d::SideSubElement(int side,int position) {
   if (position<0 || position>8 || side <0 ||side>20) {
   	PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
   }
   if(side==20) {
      return TPZGeoElSide(SubElement(position),20);
   }
   if(side<6) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>5 && side<15) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
      	int s = side-6;
         return TPZGeoElSide(SubElement(TPZCompElPr3d::RibSons[s][position]),side);
      }
   }
   if(side>14) {//faces
       if(position<0 || position>4) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
      	int s = side-15;//s = 0,1,2,3,4
         int k = 0;
         if(side==15 && position==3) k =  4;
         if(side==19 && position==3) k = -4;
         return TPZGeoElSide(SubElement(TPZCompElPr3d::FaceSons[s][position]),side+k);
      }
   }
   return TPZGeoElSide();
}
*/
/*
void TPZGeoElPr3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 20) {
      PZError << "TPZGeoElPr3d::SideSubElements called for side " << side << endl;
      return;
   }
   if(side==20) {
      sub.Resize(10);
      for(int i=0;i<8;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<6) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>5 && side<15) {//lados
      int s = side-6;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElPr3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElPr3d::SideNodes[s][1]];
      return;
   }
   if(side>14) {//faces
      int s = side-15;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElPr3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElPr3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElPr3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElPr3d::FaceSons[s][3]];
   }
}
*/
/*
TPZGeoElSide TPZGeoElPr3d::Father(int side) {

   if(!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i;
 	for(i=0;i<8;i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {
	   PZError << "TPZGeoElPr3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 6
   if(side==20) return TPZGeoElSide(fFather,side);
   //os filhos interiores não tém pai associados a seus cantos
   if(side == whichsub && side<3) return TPZGeoElSide(fFather,side);//cantos
   if((side+1) == whichsub && side>2 && side<6) return TPZGeoElSide(fFather,side);//cantos
   //lados
   if(whichsub == 0 && (side== 6 || side== 8 || side== 9)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side== 6 || side== 7 || side==10)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side== 7 || side== 8 || side==11)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side== 9 || side==12 || side==14)) return TPZGeoElSide(fFather,side);
   if(whichsub == 5 && (side==10 || side==12 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 && (side==11 || side==13 || side==14)) return TPZGeoElSide(fFather,side);
   //para os filhos 3 e 7 nenhuma aresta esta contida nalguma aresta do pai
   //faces
   if(whichsub == 0 && (side==15 || side==16 || side==18)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==15 || side==16 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==15 || side==17 || side==18)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 &&  side==19)                          return TPZGeoElSide(fFather,15);
   if(whichsub == 4 && (side==16 || side==18 || side==19)) return TPZGeoElSide(fFather,side);
   if(whichsub == 5 && (side==16 || side==17 || side==19)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 && (side==17 || side==18 || side==19)) return TPZGeoElSide(fFather,side);
   if(whichsub == 7 &&  side==15)                          return TPZGeoElSide(fFather,19);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElPr3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 5
      int s = side;
      if(side>2) s +=1;
   	if(fSubEl[s]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElPr3d::GetSubElement subelement does not contain refnode" << endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[s],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k = 0; k < refnode.NElements(); k++) {//se o subelemento k do thisside tiver o canto refnode[k]
		   for(j=0;j<6;j++) {  //este é um subelemento procurado
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   if(side ==15 || side == 19) sub[3] = SideSubElement(side,3);
   if(side == 20) {
      sub[6] = sub[5];
      sub[5] = sub[4];
      sub[4] = sub[3];
   	sub[3] = SideSubElement(side,3);
      sub[7] = SideSubElement(side,7);
   }
   return;

}
*/
/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
//transforma piramide para piramide ou pirâmide para tetraedro
/*
void TPZGeoElPr3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform &t) {
   if(!fFather || side > 20) return;
   int whichsub = -1,i;
 	for(i=0;i<8;i++)  if(fFather->SubElement(i) == this) whichsub = i;
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
   	cout << "TPZGeoElPr3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 20) {
      mult(0,0) = 0.5;//transformação entre elemento mestre do filho
      mult(1,1) = 0.5;//para o elemento mestre do pai
      mult(2,2) = 0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
         case 0:
            sum(2,0) = -.5;
            break;
         case 1:
         	sum(0,0) =  .5;
            sum(2,0) = -.5;
            break;
         case 2:
         	sum(1,0) = .5;
            sum(2,0) =-.5;
            break;
         case 3:
            mult(0,1) =  0.5;
            mult(1,1)*= -1.;
            mult(2,2)*= -1.;
            sum(1,0)  =  .5;
            sum(2,0)  = -.5;
            break;
         case 4:
            sum(2,0) =  .5;
            break;
         case 5:
         	sum(0,0) =  .5;
            sum(2,0) =  .5;
            break;
         case 6:
         	sum(1,0) =  .5;
            sum(2,0) =  .5;
            break;
         case 7:
            mult(0,1) =  0.5;
            mult(1,1)*= -1.;
            mult(2,2)*= -1.;
            sum(1,0)  =  .5;
            sum(2,0)  =  .5;
      }
   //face do filho para a face do pai
   } else if(side>14) {//15 a 19
      int s = side-15;//face do prisma
      if(whichsub==3 || whichsub==7) {
         if(side==15) s = 4;
         if(side==19) s = 0;
      }
      whichsub = -1;
      for(i=0; i<4; i++) if(fFather->SubElement(TPZCompElPr3d::FaceSons[s][i]) == this) whichsub = i;
      if(whichsub == -1) return;
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      if(side==15 || side==19) {
         switch(whichsub) {//face triangular
            case 0:
               break;
            case 1:
               sum(0,0) = 0.5;
               break;
            case 2:
               sum(1,0) = 0.5;
               break;
            case 3:
               mult(0,0) =  0.5;
               mult(0,1) =  0.5;
               mult(1,0) =  0.;
               mult(1,1) = -0.5;
               sum(0,0)  =  0.;
               sum(1,0)  =  0.5;
               break;
         }
	   } else {//face quadrilateral
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
      }
   } else if(side>5) {//e side < 15 : arestas
      whichsub = -1;
      int s = side-6;//0 a 7
      for(i=0; i<2; i++)
	      if(fFather->SubElement(TPZCompElPr3d::RibSons[s][i]) == this) whichsub = i;
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
TPZTransform TPZGeoElPr3d::SideToSideTransform(int sidefrom,int sideto) {
	return TPZShapePrism::SideToSideTransform(sidefrom,sideto);
}

TPZGeoEl *TPZGeoElPr3d::CreateBCGeoEl(int side,int bc) {
	TPZGeoEl *gel = TPZGeoPrism::CreateBCGeoEl(this,side,bc);
	return gel;
}

/*
void TPZGeoElPr3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(4,-1);
   if((face>-1 && face<6) || (face>14 && face<20)) {//cada condicao refer-se à mesma face
   	if(face>14) face = face-15;
      ids[0] = NodeIndex(TPZCompElPr3d::FaceNodes[face][0]);
      ids[1] = NodeIndex(TPZCompElPr3d::FaceNodes[face][1]);
      ids[2] = NodeIndex(TPZCompElPr3d::FaceNodes[face][2]);
		ids[3] = NodeIndex(TPZCompElPr3d::FaceNodes[face][3]);
      if(face==0 && face==4) ids.Resize(3);//face triangular
      return;
   }
 	cout << "TPZCompElPr3d::NodeFaceIds bad side , side = " << face << endl;
}
*/

// static int fatherside[8][21] = {
// /*00*/{0,6,8,9,16,18,6,15,8,9,16,18,16,20,18,15,16,20,18,20,20},
// /*01*/{6,1,7,16,10,17,6,7,15,16,10,17,16,17,20,15,16,17,20,20,20},
// /*02*/{8,7,2,18,17,11,15,7,8,18,17,11,20,17,18,15,20,17,18,20,20},
// /*03*/{18,17,16,8,7,6,20,20,20,18,17,16,15,15,15,20,20,20,20,15,20},
// /*04*/{9,16,18,3,12,14,16,20,18,9,16,18,12,19,14,20,16,20,18,19,20},
// /*05*/{16,10,17,12,4,13,16,17,20,16,10,17,12,13,19,20,16,17,20,19,20},
// /*06*/{18,17,11,14,13,5,20,17,18,18,17,11,19,13,14,20,20,17,18,19,20},
// /*07*/{14,13,12,18,17,16,19,19,19,18,17,16,20,20,20,19,20,20,20,20,20},
// };

TPZGeoElSide TPZGeoElPr3d::Father2(int side){

//	if(side<0 || side>20){
//		PZError << "TPZGeoElPr3d::Father2 called error" << endl;
//		return TPZGeoElSide();
//	}
//	if(!fFather) return TPZGeoElSide(0,0);
//	int subelindex = WhichSubel();
//	if(fatherside[subelindex][side]<0){
//		PZError << "TPZGeoElPr3d::Father2 called with index error\n";
//		return TPZGeoElSide();
//	}
//	return TPZGeoElSide(fFather,fatherside[subelindex][side]);
	if (!fFather ) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0) return TPZGeoElSide();
	return TPZGeoElSide(fFather,TPZRefPrism::FatherSide(side,son));
}

int TPZGeoElPr3d::FatherSide(int side, int son){
	return TPZRefPrism::FatherSide(side,son);
}

void TPZGeoElPr3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

	TPZRefPrism::GetSubElements(this,side,subel);
}

int TPZGeoElPr3d::NSideSubElements2(int side) {
	return TPZRefPrism::NSideSubElements(side);
}


TPZTransform TPZGeoElPr3d::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t){//Augusto:09/01/01

	if(side<0 || side>20 || !fFather){
  	PZError << "TPZGeoElPr3d::BuildTransform2 side out of range or father null\n";
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

static REAL MidSideNode[21][3] = {
/*00*/{0.,.0,-1.},/*01*/{1.,0.,-1.},/*02*/{.0,1.,-1.},/*03*/{.0,.0, 1.},
/*04*/{1.,.0, 1.},/*05*/{0.,1., 1.},/*06*/{.5,.0,-1.},/*07*/{.5,.5,-1.},
/*08*/{.0,.5,-1.},/*09*/{0.,.0, 0.},/*10*/{1.,.0, 0.},/*11*/{.0,1., 0.},
/*12*/{.5,.0, 1.},/*13*/{.5,.5, 1.},/*14*/{.0,.5, 1.},/*15*/{1./3.,1./3.,-1.},
/*16*/{.5,.0, 0.},/*17*/{.5,.5, 0.},/*18*/{0.,.5, 0.},/*19*/{1./3.,1./3., 1.},
/*20*/{1./3.,1./3.,0.} };

int TPZGeoElPr3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<8;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<21;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      TPZTransform telsd = TPZShapePrism::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform sont(son->SideDimension(sd));
      TPZTransform t = son->BuildTransform2(sd,gel,sont);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = TPZShapePrism::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(20).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao errada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(20).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(20).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}

void TPZGeoElPr3d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >7){
    PZError << "TPZGeoElPr3d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}

TPZIntPoints * TPZGeoElPr3d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoPrism::CreateSideIntegrationRule(side,order);
}

TPZTransform TPZGeoElPr3d::GetTransform(int side,int son) {
	return TPZRefPrism::GetTransform(side,son);
}

void TPZGeoElPr3d::CenterPoint(int side, TPZVec<REAL> &masscent){

  if(side < 0 || side > NSides()-1){
    PZError << "TPZGeoElPr3d::CenterPoint error side = " << side << endl;
    return;
  }
  TPZShapePrism::CenterPoint(side,masscent);
}

REAL TPZGeoElPr3d::RefElVolume(){
  return TPZShapePrism::RefElVolume();
}
