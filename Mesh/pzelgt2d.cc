//$Id: pzelgt2d.cc,v 1.6 2003-11-05 16:02:21 tiago Exp $

// -*- c++ -*-
//METHODS DEFINITION FOR CLASS ELEMT2D

#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelgq2d.h"
#include "pzelg1d.h"
#include "pzshapetriang.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzgnode.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzelgpoint.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"

static TPZCompEl *CreateEl(TPZGeoElT2d *gel,TPZCompMesh &mesh,int &index) {
  return new TPZCompElT2d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElT2d::fp)(TPZGeoElT2d *,TPZCompMesh &,int &) = CreateEl;

/**The vector with nodes has first the three corners nodes. The last three nodes are for middle nodes
   on the sides*/
TPZGeoElT2d::TPZGeoElT2d(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh)
  : TPZGeoEl(id,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod != 3) {
    PZError << "TPZGeoElT2d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<3;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElT2d::TPZGeoElT2d() {
  int i;
  for(i=0;i<3;i++) fNodeIndexes[i] = -1;
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElT2d::TPZGeoElT2d(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {


  int i,nnod = nodeindexes.NElements();
  if(nnod != 3) {
    PZError << "TPZGeoElT2d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<3;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElT2d::~TPZGeoElT2d() { }

void TPZGeoElT2d::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
	TPZGeoTriangle::Shape(param,phi,dphi);
}

TPZGeoElT2d *TPZGeoElT2d::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElT2d(nodeindexes,matind,mesh);
}

int TPZGeoElT2d::NodeIndex(int node) {
  if(node<0 || node>2) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElT2d::NSideNodes(int side) {
  if(side<0 || side>6) {
    PZError << "TPZGeoElT2d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<3) return 1;
  if(side<6) return 2;
  return 3;
}

int TPZGeoElT2d::SideNodeIndex(int side,int node) {
	int loc = TPZShapeTriang::SideNodeLocId(side,node);
	if(loc>=0) return fNodeIndexes[loc];
	return -1;
}

int TPZGeoElT2d::SideNodeLocIndex(int side,int node) {
	return TPZShapeTriang::SideNodeLocId(side,node);
}

void TPZGeoElT2d::MidSideNodeIndex(int side,int &index) {
  index = -1;
  if(side<0 || side>5) {
    PZError << "TPZGeoElT2d::MidSideNodeIndex. Bad parameter side = " << side << endl;
    return;
  }
  TPZRefTriangle::MidSideNodeIndex(this,side,index);
}

void TPZGeoElT2d::NewMidSideNode(int side,int &index) {
	TPZRefTriangle::NewMidSideNode(this,side,index);
}

/**Determine the coordinate of the center of the element*/
/*
int TPZGeoElT2d::CenterIndex() {
  TPZVec<REAL> coord(3);
  TPZVec<REAL> param(2,1./3.);
  X(param,coord);  //It determine the centroid of the triangular element.
  int indexcenter = Mesh()->NodeVec().AllocateNewElement();
  Mesh()->NodeVec()[indexcenter].Initialize(coord,*Mesh());
  return indexcenter;
  // We can not to take the center type of the master cell because of characteristic
  // point is not preserved for the transformation.
}
*/
int TPZGeoElT2d::SideDimension(int side) {
	return TPZShapeTriang::SideDimension(side);
}

/*
TPZGeoElSide TPZGeoElT2d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 ou 2
//o lado 6 tem dim. 2 e neste caso targetdimension = 3 logo : 0 <= side <= 5
  if( (side<0 || side>5) || (targetdimension != 1 && targetdimension != 2) ) {
    PZError << "TPZGeoElT2d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
//side = 0,1,2,3,4,5
  switch(targetdimension) {
	  case 1:
       if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,3);
       	 if(side==2) return TPZGeoElSide(this,5);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,3);
       	 if(side==2) return TPZGeoElSide(this,4);
       } else if(this == father->SubElement(2)) {
       	 if(side==0) return TPZGeoElSide(this,5);
       	 if(side==1) return TPZGeoElSide(this,4);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(father->SubElement(1),4);
       	 if(side==1) return TPZGeoElSide(father->SubElement(0),5);
       	 if(side==2) return TPZGeoElSide(father->SubElement(0),3);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
       	return TPZGeoElSide(this,6);
  }//switch
  return TPZGeoElSide();
}
*/
void TPZGeoElT2d::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> high;
	TPZShapeTriang::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}
}

void TPZGeoElT2d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapeTriang::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapeTriang::SideConnectLocId(side,is));

}

/*
void TPZGeoElT2d::SideMasterCo(int side,TPZFMatrix &coord) {
// CEDRIC VERIFICAR , comparar com 1D : verificado!!!
  if(side<0 || side>6) {
    PZError << "TPZGeoElT2d::SideMasterCo. Bad parameter side.\n";
    return;
  }
  int row = coord.Rows();//2x3
  if(side==6) {
    coord.Redim(row,3);          //side 6 tem dimensão 2
    coord(0,1) = coord(1,2) = 1.;// [ {0,1,0},{0,0,1} ]
    return;                      //ou (0,0),(1,0),(0,1) nós 0,1,2
  }
  if(side<3) {         //side<3 tem dimensão 0
    coord.Redim(row,1);//2x1
    coord.Zero();      //side = 0 => (0,0)
    if(side==1) coord(0,0)=1.;//(1,0)
    else if(side==2) coord(1,0)=1.;//(0,1)
    return;
  }
  coord.Redim(row,2);  //side = 3,4,5 tem dimensão 1
  coord.Zero();        //2x2  [{0,0},{0,0}]
  if(side==3) coord(0,1) = 1.;//(0,0),(1,0)
  else if(side==5) coord(1,0) = 1.;//(0,1),(0,0)
  else if (side == 4) coord(0,0) = coord(1,1) = 1.;//(1,0),(0,1)
}
*/
void TPZGeoElT2d::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
	TPZFMatrix xco(3,TPZShapeTriang::NNodes);
	int i,j;
	for(i=0; i<TPZShapeTriang::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoTriangle::Jacobian(xco,param,jacobian,axes,detjac,jacinv);
}

void TPZGeoElT2d::X(TPZVec<REAL> &par,TPZVec<REAL> &result){
	TPZFMatrix xco(3,TPZGeoTriangle::NNodes);
	int i,j;
	for(i=0; i<TPZShapeTriang::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoTriangle::X(xco,par,result);
}

/*
void TPZGeoElT2d::NormalVector(int side,TPZVec<REAL> &loc,
			     TPZVec<REAL> &normal,TPZFMatrix &axes,TPZFMatrix &jacside) {
  if(side < 0 || side > 6) {
    PZError << "TPZGeoElT2d.NormalVector invalid side : "
      " side = " << side << "\n";
    return;
  }

  TPZVec<REAL> t(3,0.);
  int i,id,ic,j=(side+3)%3;
  REAL detjac;
  TPZFMatrix jacinv(0);
  if(side==6 || side<3) {
    REAL norm = 0.;
    TPZVec<REAL> t1(3,0.);
    for(i=0;i<3;i++) {
      t[i]= NodePtr(j)->Coord(i);
      t1[i]= t[i] - NodePtr((j+1)%3)->Coord(i);
      t[i]-=NodePtr((j+2)%3)->Coord(i);
    }
    if(side==6) {
    	normal[0]=t[1]*t1[2]-t1[1]*t[2];
      normal[1]=t[2]*t1[0]-t1[2]*t[0];
      normal[2]=t[0]*t1[1]-t[1]*t1[0];
      jacside.Redim(2,2);
      Jacobian(loc,jacside,axes,detjac,jacinv);
    }
    else {
      for(i=0;i<3;i++) normal[i]=t[i]+t1[i];
      jacside.Redim(0,0);
      axes.Zero();
      for(i=0;i<3;i++) axes(i,i)=1.;
    }
    for(i=0;i<3;i++) norm += normal[i]*normal[i];
    norm = sqrt(norm);
    for(i=0;i<3;i++) normal[i]/=norm;
    return;
  }
  REAL spacephi[3],spacedphi[6];
  TPZFMatrix phi(3,1,spacephi,3);
  TPZFMatrix dphi(2,3,spacedphi,6);
  Shape(loc,phi,dphi);

  TPZGeoNode* np;
  REAL ider[2] = {0.,0.}, sq2 = sqrt(2.)/2.;
  switch(j) {
  case 0: ider[0] = 1.;break;
  case 1: ider[0] = -sq2; ider[1] = sq2; break;
  case 2: ider[1] = -1.; break;
  }
  for(i=0;i<3;i++) {
    np = NodePtr(i);
    for(ic=0;ic<3;ic++)
      for(id=0;id<2;id++)
	     t[ic] += ider[id] * (np->Coord(ic)) * dphi(id,i);
  }

  REAL jac1dvar,tnorm;
  tnorm = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
  if(j != 1) jac1dvar = 0.5*tnorm;
  else jac1dvar = sq2*tnorm;
  jacside(0,0) = jac1dvar;

  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.),V1til(3,0.),V1Ttil(3,0.);
  REAL V1Norm=0.,V2Norm=0.,V1V2=0.,V2tilNorm=0.,V1tilNorm =0.,tV1 = 0.,V1TtilNorm = 0.;
  for(i=0;i<3;i++) {
    np = NodePtr(i);
    for(ic=0;ic<3;ic++) {
      V1[ic] += np->Coord(ic)*dphi(0,i);
      V2[ic] += np->Coord(ic)*dphi(1,i);
    }
  }
  for(ic=0;ic<3;ic++) {
    V1Norm += V1[ic]*V1[ic];
    V2Norm += V2[ic]*V2[ic];
    V1V2 += V1[ic]*V2[ic];
    tV1 += t[ic]*V1[ic];
  }
  V1Norm = sqrt(V1Norm);
  V2Norm = sqrt(V2Norm);
  for(ic=0;ic<3;ic++) {
    V1[ic] /= V1Norm;
    V2[ic] /= V2Norm;
    //Jorge 14/10/99
    t[ic] /= tnorm;
    V2til[ic] = V2[ic] - V1V2*V1[ic]/V1Norm/V2Norm;
    V1til[ic] = V1[ic] - V1V2*V2[ic]/V1Norm/V2Norm;
    //Jorge 14/10/99
    V1Ttil[ic] = V1[ic] - tV1*t[ic]/V1Norm/tnorm;
    V2tilNorm += V2til[ic]*V2til[ic];
    V1tilNorm += V1til[ic]*V1til[ic];
    V1TtilNorm += V1Ttil[ic]*V1Ttil[ic];
  }
  V2tilNorm = sqrt(V2tilNorm);
  V1tilNorm = sqrt(V1tilNorm);
  V1TtilNorm = sqrt(V1TtilNorm);
  for(ic=0;ic<3;ic++) {
    axes(0,ic) = V1[ic];
    axes(1,ic) = V2til[ic]/V2tilNorm;
  }
  switch(j) {
  case 0:
    normal[0] = -V2til[0]/V2tilNorm;
    normal[1] = -V2til[1]/V2tilNorm;
    normal[2] = -V2til[2]/V2tilNorm;
    break;
  case 1:
    normal[0] = V1Ttil[0]/V1TtilNorm;
    normal[1] = V1Ttil[1]/V1TtilNorm;
    normal[2] = V1Ttil[2]/V1TtilNorm;
    break;
  case 2:
    normal[0] = -V1til[0]/V1tilNorm;
    normal[1] = -V1til[1]/V1tilNorm;
    normal[2] = -V1til[2]/V1tilNorm;
    break;
  }
  axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
  axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
  axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
}
*/
/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElT2d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
	TPZRefTriangle::Divide(this,SubElVec);
}

int TPZGeoElT2d::NSubElements() {
	return TPZRefTriangle::NSubEl;
}

//int TPZGeoElT2d::NSideSubElements(int side) {
//	return TPZRefTriangle::NSideSubElements(side);
//}

/*
void TPZGeoElT2d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
  if(!fSubEl[0]) {
    sub.Resize(0);
    return;
  }
  if(side < 0 || side > 6) {
    PZError << "TPZGeoElT2d::SideSubElements called for side " << side << endl;
    return;
  }
  if(side==6) {//6
    for(int i=0;i<4;i++) sub[i] = fSubEl[i];
    return;
  }
  //0,1,2
  if(side<3) {
    sub[0]=fSubEl[side];
    return;
  }
  //3,4,5
  side-=3;
  sub.Resize(2);
  sub[0] = fSubEl[side];
  sub[1] = fSubEl[(side+1)%3];
}
*/
/*
TPZGeoElSide TPZGeoElT2d::SideSubElement(int side,int position) {
   if (position<0 ||position>3 || side <0 ||side>6) {
   	PZError << "TPZGeoElT2d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
   }
   if(side==6) return TPZGeoElSide(SubElement(position),6);
   if(side<3) {
      if(position!=0) {
         PZError << "TPZGeoElT2d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(position==0) return TPZGeoElSide(SubElement(side-3),side);
   else if(position == 1) return TPZGeoElSide(SubElement((side-2)%3),side);
   PZError << "TPZGeoElT2d::SideSubElement called with position " << position << " side " << side << endl;
   return TPZGeoElSide();
}
*/
/*
void TPZGeoElT2d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {
  int nsub = NSideSubElements(side);
  sub.Resize(nsub);
  if(!nsub) return;
  if(side<3 && refnode[0]==NodeIndex(side)) {
    sub[0]=TPZGeoElSide(fSubEl[side],side);
    return;
  }
//  if(nsub==1) {  // then side=0 or =1 or =2
//    sub[0]=TPZGeoElSide(fSubEl[side],side);
//    return;
//  }

  int i,j,k,nref=refnode.NElements();
  for(i=0;i<nsub;i++) {
    TPZGeoElSide sidesub = SideSubElement(side,i);
    TPZGeoEl *subel = sidesub.Element();
    for(k=0; k<nref; k++) {
      for(j=0;j<3;j++) {
         if(subel->NodeIndex(j)==refnode[k]) {
            sub[k] = SideSubElement(side,i);
         }
      }
    }
  }
  if(side==6) sub[3] = SideSubElement(side,3);
  return;
}
*/
/*
TPZGeoElSide TPZGeoElT2d::Father(int side) {

   if(!fFather) return TPZGeoElSide();

   int whichsub = -1;
   int i;
   for(i=0; i<4; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {//equivale a is = 4  ou is > 3
	   PZError << "TPZGeoElT2d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 4
   if(whichsub == 3 && side != 6) return TPZGeoElSide();
   if(whichsub == side || side==6) return TPZGeoElSide(fFather,side);//side = 0,1,2
	if(whichsub != side && side<3) return TPZGeoElSide();
   //side = 3,4,5
 	if(whichsub == 0 && side!=4) return TPZGeoElSide(fFather,side);//ou side==3 || side==5
 	if(whichsub == 1 && side!=5) return TPZGeoElSide(fFather,side);//ou side==3 || side==4
 	if(whichsub == 2 && side!=3) return TPZGeoElSide(fFather,side);//ou side==4 || side==5
 	//if(wichsub == 3) return TPZGeoElSide();//é feito pelo seguinte caso

//   PZError << "TPZGeoElT2d::Father. fFather isn't father of this element along the given side.\n";
   return TPZGeoElSide();//inclui os outros casos
}
*/
/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
/*
void TPZGeoElT2d::BuildTransform(int side,TPZGeoEl *father,TPZTransform &t) {
  if(this == father) return;
   if(!fFather) {
     cout << "TPZGeoElT2d::BuildTransform called for inconsistent parameters\n";
     return;
   }
   int whichsub = -1;
   int i;
   for(i=0; i<4; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) return;
   int dim = SideDimension(side);
   TPZTransform tloc(dim);
   REAL store[6];
   TPZFMatrix mult(dim,dim,store,4);
   TPZFMatrix sum(dim,1,store+4,2);
   mult.Zero();
   sum.Zero();

   TPZGeoEl *locfather;
   if(father == this) {
	   mult(0,0) = mult(1,1) = 1.;
      locfather = fFather;
   } else {
	   locfather = fFather;
   }
   if(!locfather) {
   	cout << "TPZGeoElT2d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 6) {
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      //  if(typediv==1) {
      switch(whichsub) {
         case 0:
            break;
         case 1:
            sum(0,0) = 0.5;
            break;
         case 2:
            sum(1,0) = 0.5;
            break;
         case 3:
            mult(0,0) = -0.5;
            mult(0,1) = 0.;
            mult(1,0) = 0.;
            mult(1,1) = -0.5;
            sum(0,0) = 0.5;
            sum(1,0) = 0.5;
            break;
      }
   } else if(side <3) {
   	return;
   } else {
   	mult(0,0) = 0.5;
      if(whichsub == side-3) sum(0,0) = -0.5;
      else sum(0,0) = 0.5;
   }

   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform(side,father,t);
}
*/
TPZTransform TPZGeoElT2d::SideToSideTransform(int sidefrom,int sideto) {
	return TPZShapeTriang::SideToSideTransform(sidefrom,sideto);
}

TPZGeoEl *TPZGeoElT2d::CreateBCGeoEl(int side,int bc) {
	TPZGeoEl *geo = TPZGeoTriangle::CreateBCGeoEl(this,side,bc);
	return geo;
}

TPZGeoElSide TPZGeoElT2d::Father2(int side){//Augusto:09/01/01
	if (!fFather ) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0 || !fFather) return TPZGeoElSide();
	return TPZGeoElSide(fFather,TPZRefTriangle::FatherSide(side,son));
}

int TPZGeoElT2d::FatherSide(int side, int son){
	return TPZRefTriangle::FatherSide(side,son);
}

void TPZGeoElT2d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){//Augusto:09/01/01

	TPZRefTriangle::GetSubElements(this,side,subel);
}

int TPZGeoElT2d::NSideSubElements2(int side) {
	return TPZRefTriangle::NSideSubElements(side);
}


TPZTransform TPZGeoElT2d::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t){//Augusto:09/01/01


	if(side<0 || side>TPZShapeTriang::NSides || !fFather){
  	PZError << "TPZGeoElT2d::BuildTransform2 side out of range or father null\n";
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

static REAL MidSideNode[7][3] = {
/*00*/{.0,0.},/*01*/{1.0,.0},/*02*/{0.,1.0},
/*03*/{.5,0.},/*04*/{0.5,.5},/*05*/{0.,0.5},
/*06*/{ 1./3.,1./3.} };

int TPZGeoElT2d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  pss[2] = 0.;//2d
  pfs[2] = 0.;
  pf[2] = 0.;
  for(sn=0;sn<4;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<7;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];// = 0
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      TPZTransform telsd = TPZShapeTriang::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform sont(son->SideDimension(sd));
      TPZTransform t = son->BuildTransform2(sd,gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = TPZShapeTriang::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(6).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(6).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(6).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}

/** Compute the measure of the geometrical element - Jorge 17/7/99*/
/*
REAL TPZGeoElT2d::Mesure(int dim) {
  if(dim!=2) return 0.;
  if(fMesure == 0.) {
    TPZGeoNode &nod1 = Mesh()->NodeVec()[fNodeIndexes[0]];
    REAL xx, yy, x0 = nod1.Coord(0), y0 = nod1.Coord(1);
    TPZGeoNode &nod2 = Mesh()->NodeVec()[fNodeIndexes[1]];
    xx = x0 - nod2.Coord(0);
    yy = nod2.Coord(1) - y0;
    TPZGeoNode &nod3 = Mesh()->NodeVec()[fNodeIndexes[2]];
    xx *= (nod3.Coord(1) - y0);
    yy *= (nod3.Coord(0) - x0);
    fMesure = 0.5 * fabs(xx + yy);
  }
  return fMesure;
}
*/
//void TPZGeoElT2d::Center(TPZVec<REAL> &center) {
//  center[0] = center[1] = 1./3.;
//}

void TPZGeoElT2d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >3){
    PZError << "TPZGeoElT2d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}

TPZIntPoints * TPZGeoElT2d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoTriangle::CreateSideIntegrationRule(side,order);
}

TPZTransform TPZGeoElT2d::GetTransform(int side,int son) {
	return TPZRefTriangle::GetTransform(side,son);
}

void TPZGeoElT2d::CenterPoint(int side, TPZVec<REAL> &masscent){

  TPZShapeTriang::CenterPoint(side,masscent);
}

REAL TPZGeoElT2d::RefElVolume(){
  return TPZShapeTriang::RefElVolume();
}
