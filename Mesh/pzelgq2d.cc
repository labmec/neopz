//$Id: pzelgq2d.cc,v 1.4 2003-11-05 16:02:21 tiago Exp $

// -*- c++ -*-
//METHODS DEFINITION FOR CLASS ELEMQ2D

#include "pzelgq2d.h"
#include "pzelg1d.h"
#include "pzelcq2d.h"
#include "pzelgpoint.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzshapequad.h"
#include "pzrefquad.h"
#include "pzgeoquad.h"

static TPZCompEl *CreateEl(TPZGeoElQ2d *gel,TPZCompMesh &mesh,int &index) {
  return new TPZCompElQ2d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElQ2d::fp)(TPZGeoElQ2d *,TPZCompMesh &,int &) = CreateEl;

TPZGeoElQ2d::TPZGeoElQ2d(int id,TPZVec<int> &nodeindexes,int matid,TPZGeoMesh &mesh):
  TPZGeoEl(id,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElQ2d::Constuctor, number of nodes : " << nnod << endl;
    if(nnod != 9) return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElQ2d::TPZGeoElQ2d(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElQ2d::Constuctor, number of nodes : " << nnod << endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElQ2d::~TPZGeoElQ2d() {}

void TPZGeoElQ2d::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {

   REAL x=param[0], y=param[1];

   phi(0,0) = .25*(1.-x)*(1.-y);
   phi(1,0) = .25*(1.+x)*(1.-y);
   phi(2,0) = .25*(1.+x)*(1.+y);
   phi(3,0) = .25*(1.-x)*(1.+y);

   dphi(0,0) = .25*(y-1.);
   dphi(1,0) = .25*(x-1.);

   dphi(0,1) = .25*(1.-y);
   dphi(1,1) =-.25*(1.+x);

   dphi(0,2) = .25*(1.+y);
   dphi(1,2) = .25*(1.+x);

   dphi(0,3) =-.25*(1.+y);
   dphi(1,3) = .25*(1.-x);


/*   REAL spacephi1[3],spacephi2[3],spacedphi1[3],spacedphi2[3];
   TPZFMatrix phi1(2,1,spacephi1,3);
   TPZFMatrix phi2(2,1,spacephi2,3);
   TPZFMatrix dphi1(1,2,spacedphi1,3);
   TPZFMatrix dphi2(1,2,spacedphi2,3);
   REAL x=param[0], y=param[1];

   TPZGeoEl::Shape1d(x,2,phi1,dphi1);
   TPZGeoEl::Shape1d(y,2,phi2,dphi2);

   phi(0,0) = phi1(0,0)*phi2(0,0);
   dphi(0,0) = dphi1(0,0)*phi2(0,0);
   dphi(1,0) = phi1(0,0)*dphi2(0,0);
   phi(1,0) = phi1(1,0)*phi2(0,0);
   dphi(0,1) = dphi1(0,1)*phi2(0,0);
   dphi(1,1) = phi1(1,0)*dphi2(0,0);
   phi(2,0) = phi1(1,0)*phi2(1,0);
   dphi(0,2) = dphi1(0,1)*phi2(1,0);
   dphi(1,2) = phi1(1,0)*dphi2(0,1);
   phi(3,0) = phi1(0,0)*phi2(1,0);
   dphi(0,3) = dphi1(0,0)*phi2(1,0);
   dphi(1,3) = phi1(0,0)*dphi2(0,1);*/
}

TPZGeoElQ2d *TPZGeoElQ2d::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElQ2d(nodeindexes,matind,mesh);
}

int TPZGeoElQ2d::NNodes() {
	return TPZGeoQuad::NNodes;
}
int TPZGeoElQ2d::NodeIndex(int node) {
  if(node<0 || node>3) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElQ2d::NSideNodes(int side) {
	return TPZShapeQuad::NSideNodes(side);
}

int TPZGeoElQ2d::SideNodeIndex(int side,int node) {
  if(side<0 || side>8) {
    PZError << "TPZGeoElQ2d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return fNodeIndexes[TPZShapeQuad::SideNodeLocId(side,node)];
}

int TPZGeoElQ2d::SideNodeLocIndex(int side,int node) {
	return TPZShapeQuad::SideNodeLocId(side,node);
}

void TPZGeoElQ2d::MidSideNodeIndex(int side,int &index) {
	TPZRefQuad::MidSideNodeIndex(this,side,index);
}

void TPZGeoElQ2d::NewMidSideNode(int side,int &index) {
	TPZRefQuad::NewMidSideNode(this,side,index);
}

/**Determine the coordinate of the center of the element*/
/*
int TPZGeoElQ2d::CenterIndex() {
  TPZGeoElSide nghsd = Neighbour(8);
  while(nghsd.Exists() && !nghsd.HasSubElement() && nghsd.Element() != this) nghsd = nghsd.Neighbour();
// Philippe 15/4/99
  if(nghsd.Exists() && nghsd.Element()->HasSubElement()) {
	  int index = -1;
	  nghsd.Element()->MidSideNodeIndex(nghsd.Side(),index);
	  return index;
  }
  TPZVec<REAL> coord(3);
  TPZVec<REAL> param(2,0);
  X(param,coord);//It determine the centroid of the quadrilateral element.
  int indexcenter = Mesh()->NodeVec().AllocateNewElement();
  Mesh()->NodeVec()[indexcenter].Initialize(coord,*Mesh());
  return indexcenter;
  // We can not to take the center type of the master cell because of characteristic
  // point is not preserved for the transformation.
}
*/
int TPZGeoElQ2d::SideDimension(int side) {
	return TPZShapeQuad::SideDimension(side);
}
/*
TPZGeoElSide TPZGeoElQ2d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 ou 2
//se side = 0,1,2,3 targetdimension deve ser 1
//se side = 4,5,6,7 targetdimension deve ser 2
  if( (side<0 || side>7) || (targetdimension != 1 && targetdimension != 2) ) {
    PZError << "TPZGeoElQ2d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension << endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
//side = 0,1,2,3,4,5,6,7
  switch(targetdimension) {//=1,2
	  case 1:
       if       (this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,4);
       	 if(side==3) return TPZGeoElSide(this,7);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,4);
       	 if(side==2) return TPZGeoElSide(this,5);
       } else if(this == father->SubElement(2)) {
       	 if(side==1) return TPZGeoElSide(this,5);
       	 if(side==3) return TPZGeoElSide(this,6);
       } else if(this == father->SubElement(3)) {
       	 if(side==2) return TPZGeoElSide(this,6);
       	 if(side==0) return TPZGeoElSide(this,7);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
       	return TPZGeoElSide(this,8);
  }//switch
  return TPZGeoElSide();
}
*/
void TPZGeoElQ2d::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> high;
	TPZShapeQuad::HigherDimensionSides(side,high);
	int cap = high.NElements(),s;
	for(s=0; s<cap; s++) {
		if(SideDimension(high[s]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,high[s]));
		}
	}

}


void TPZGeoElQ2d::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	int nsidecon = TPZShapeQuad::NSideConnects(side);
	int is;
	for(is=0; is<nsidecon-1; is++) smallsides.Push(TPZShapeQuad::SideConnectLocId(side,is));
}

/*
void TPZGeoElQ2d::SideMasterCo(int side,TPZFMatrix &coord) {
  if(side<0 || side>8) {
    PZError << "TPZGeoElQ2d::SideMasterCo. Bad parameter side.\n";
    return;
  }
  int row = coord.Rows();
  if(side==8) {//side 8 tem dimensão 2
    coord.Redim(row,4);//2x4
    coord(0,1) = coord(1,2) = 1.;// [ {0,1,0,-1},{-1,0,1,0} ]
    coord(0,3) = coord(1,0) =-1.;//ou (0,-1),(1,0),(0,1),(-1,0) nós 0,1,2,3
    return;
  }
  if(side<3) {         //side=0,1,2,3 tem dimensão 0
    coord.Redim(row,1);//2x1 um unico nó
    if     (side==0) {coord(0,0)=-1.;coord(1,0)=-1.;}//(-1,-1) : 0
    else if(side==1) coord(0,0)= 1.;//(1,-1) : 1
    else if(side==2) coord(1,0)= 1.;//(1,1) : 2
    else if(side==3) coord(0,0)=-1.;//(-1,1) : 3
    return;
  }
  coord.Redim(row,2);  //side = 4,5,6,7 tem dimensão 1
  coord.Zero();        //2x2  [{0,0},{0,0}]
  if     (side==4) {coord(0,0)=-1.;coord(0,1)= 1.;coord(1,0)=-1.;coord(1,1)=-1.;}//(-1,-1),(1,-1) : 0,1
  else if(side==5) {coord(0,0)= 1.;coord(0,1)= 1.;coord(1,1)= 1.;}//(1,-1),(1,1) : 1,2
  else if(side==6) {coord(1,0)= 1.;coord(0,1)=-1.;}//(1,1),(-1,1) : 2,3
  else if(side==7) {coord(0,0)=-1.;coord(1,1)=-1.;}//(1,1),(-1,1) : 3,0
}
*/
void TPZGeoElQ2d::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
	TPZFMatrix xco(3,TPZGeoQuad::NNodes);
	int i,j;
	for(i=0; i<TPZShapeQuad::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoQuad::Jacobian(xco,param,jacobian,axes,detjac,jacinv);
}

void TPZGeoElQ2d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
	TPZFMatrix xco(3,TPZGeoQuad::NNodes);
	int i,j;
	for(i=0; i<TPZShapeQuad::NNodes; i++) {
		for(j=0; j<3; j++) {
			xco(j,i) = NodePtr(i)->Coord(j);
		}
	}
	TPZGeoQuad::X(xco,loc,result);
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/
/*
void TPZGeoElQ2d::NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,
			     TPZFMatrix &axes,TPZFMatrix &jac1d) {
  int nnodes = NNodes();
#ifdef DEBUG
  if (nnodes != 4) {
    PZError << "TPZGeoElQ2d.NormalVector, only implemented for"
      " 4 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 2 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1.) {
    PZError << "TPZGeoElQ2d.NormalVector, fl out of range : "
      " point.NElements() = " << param.NElements() <<
      "\npoint[0] = " << param[0] << " point[1] = " << param[1] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "elgq2d.NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side >= 4) {
    PZError << "TPZGeoElQ2d.jacobian invalid side : "
      " side = " << side << "\n";
    return;
  }
#endif

  REAL spacephi[4],spacedphi[8];
//  TPZFMatrix phi(4,1,spacedphi,4);
// Philippe 31;3;99
  TPZFMatrix phi(4,1,spacephi,4);
  TPZFMatrix dphi(2,4,spacedphi,8);
  Shape(param,phi,dphi);
  TPZGeoNode *np;
  TPZVec<REAL> t(3,0.);
  int i,j,ider = 0;
  if(side==1 || side==3) ider = 1;

  for(i=0;i<nnodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++)
      t[j] += np->Coord(j)*dphi(ider,i);
  }

  //      note that ||t|| != 1 , ||t|| = |J|

  jac1d(0,0) = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );

  // consistent axes computation Philippe 17/4/97
  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.),V1til(3,0.);
  REAL V1Norm=0.,V2Norm=0.,V1V2=0.,V2tilNorm=0.,V1tilNorm =0.;
  for(i=0;i<nnodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++) {
      V1[j] += np->Coord(j)*dphi(0,i);
      V2[j] += np->Coord(j)*dphi(1,i);
    }
  }
  for(j=0;j<3;j++) {
    V1Norm += V1[j]*V1[j];
    V2Norm += V2[j]*V2[j];
    V1V2 += V1[j]*V2[j];
  }
  V1Norm = sqrt(V1Norm);
  V2Norm = sqrt(V2Norm);
  for(j=0;j<3;j++) {
    V1[j] /= V1Norm;
    V2[j] /= V2Norm;
    V2til[j] = V2[j] - V1V2*V1[j]/V1Norm/V2Norm;
    V1til[j] = V1[j] - V1V2*V2[j]/V1Norm/V2Norm;
    V2tilNorm += V2til[j]*V2til[j];
    V1tilNorm += V1til[j]*V1til[j];
  }
  V2tilNorm = sqrt(V2tilNorm);
  V1tilNorm = sqrt(V1tilNorm);
  for(j=0;j<3;j++) {
    axes(0,j) = V1[j];
    axes(1,j) = V2til[j]/V2tilNorm;
  }
  switch(side) {
  case 0:
  case 2:
    for(i=0;i<3;i++)
      normal[i] = V2til[i]/V2tilNorm;
    break;
  case 1:
  case 3:
    for(i=0;i<3;i++)
      normal[i] = V1til[i]/V1tilNorm;
    break;
  }
  if(side == 0 || side == 3) for(i=0;i<3;i++) normal[i] *= -1.;

  axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
  axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
  axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
  return;
}
*/
/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElQ2d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
   //If exist fSubCell return this sons
	TPZRefQuad::Divide(this,SubElVec);
}

/*
TPZGeoElSide TPZGeoElQ2d::SideSubElement(int side,int position) {
   if (position<0 ||position>3 || side <0 ||side>8) {
   	PZError << "TPZGeoElQ2d::SideSubElement called with position " << position << " side " << side << endl;
      return TPZGeoElSide();
   }
   if(side==8) return TPZGeoElSide(SubElement(position),8);
   if(side<4) {
      if(position!=0) {
         PZError << "TPZGeoElQ2d::SideSubElement called with position " << position << " side " << side << endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   //side = 4,5,6,7
   if(position==0) return TPZGeoElSide(SubElement(side-4),side);
   else return TPZGeoElSide(SubElement((side-3)%4),side);

}
*/
/*
void TPZGeoElQ2d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
  if(!fSubEl[0]) {
    sub.Resize(0);
    return;
  }
  if(side < 0 || side > 8) {
    PZError << "TPZGeoElQ2d::SideSubElements called for side " << side << endl;
    return;
  }
  if(side==8) {//8
    for(int i=0;i<4;i++) sub[i] = fSubEl[i];
    return;
  }
  //0,1,2,3
  if(side<4) {
    sub[0]=fSubEl[side];
    return;
  }
  //4,5,6,7
  side-=4;
  sub.Resize(2);
  sub[0] = fSubEl[side];      //0,1,2,3
  sub[1] = fSubEl[(side+1)%4];//1,2,3,0
}
*/
/*
TPZGeoElSide TPZGeoElQ2d::Father(int side) {

   if(!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i;
   for(i=0; i<4; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {//equivale a is = 4  ou is > 3
	   PZError << "TPZGeoElQ2d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 4
   if(whichsub == side || side==8) return TPZGeoElSide(fFather,side);//side = 0,1,2,3
   //side = 4,5,6,7
 	if(whichsub == 0 && (side==4 || side==7)) return TPZGeoElSide(fFather,side);//ou
 	if(whichsub == 1 && (side==4 || side==5)) return TPZGeoElSide(fFather,side);//ou side==3 || side==4
 	if(whichsub == 2 && (side==5 || side==6)) return TPZGeoElSide(fFather,side);//ou side==4 || side==5
 	if(whichsub == 3 && (side==6 || side==7)) return TPZGeoElSide(fFather,side);//ou side==4 || side==5
 	//if(wichsub == 3) return TPZGeoElSide();//é feito pelo seguinte caso

//   PZError << "TPZGeoElQ2d::Father. fFather isn't father of this element along the given side.\n";
   return TPZGeoElSide();//inclui os outros casos

}
*/
/*
void TPZGeoElQ2d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {
// REVISAR
   int nsub = NSideSubElements(side);
   sub.Resize(nsub);
   if(!nsub) return;

   int i,j,k;
   if(nsub==1) {//side = 0,1,2,3
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElQ2d::GetSubElement subelement does not contain refnode" << endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //nsub = 4
   //if(nsub == 4) nsub = 3; se nsub = 4 => side = 8
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k=0; k<nsub; k++) {
		   for(j=0;j<4;j++) {
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   //if(side == 8) sub[3] = TPZGeoElSide(fSubEl[3],8);
   return;
}
*/
/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
/*
void TPZGeoElQ2d::BuildTransform(int side,TPZGeoEl *father,TPZTransform &t) {
  if(father == this) return;
   if(!fFather) {
     cout << "TPZGeoElQ2d::BuildTransform called for inconsistent parameters\n";
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
	   locfather = fFather;//pai do lemento atual
   }
   if(!locfather) {
   	PZError << "TPZGeoElT2d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 8) {//face para face
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
      }
   } else if(side > 3 && side < 8) {
   	mult(0,0) = 0.5;
      if(whichsub == side-4) sum(0,0) = -0.5;
      else sum(0,0) = 0.5;
   }

   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform(side,father,t);
}
*/
TPZTransform TPZGeoElQ2d::SideToSideTransform(int sidefrom,int sideto) {
	return TPZShapeQuad::SideToSideTransform(sidefrom,sideto);
}

TPZGeoEl *TPZGeoElQ2d::CreateBCGeoEl(int side,int bc) {
	TPZGeoEl *gel = TPZGeoQuad::CreateBCGeoEl(this,side,bc);
	return gel;
}

TPZGeoElSide TPZGeoElQ2d::Father2(int side){//Augusto:09/01/01
	if(!fFather) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0) return TPZGeoElSide();
	return TPZGeoElSide(fFather,TPZRefQuad::FatherSide(side,son));

}

int TPZGeoElQ2d::FatherSide(int side, int son){
	return TPZRefQuad::FatherSide(side,son);
}

void TPZGeoElQ2d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){//Augusto:09/01/01
	TPZRefQuad::GetSubElements(this,side,subel);
}

int TPZGeoElQ2d::NSideSubElements2(int side) {
	return TPZRefQuad::NSideSubElements(side);
}

TPZTransform TPZGeoElQ2d::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t){//Augusto:09/01/01
	
	if(side<0 || side>(TPZShapeQuad::NSides-1) || !fFather){
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

static REAL MidSideNode[9][3] = {
/*00*/{-1.,-1.},/*01*/{ 1.,-1.},/*02*/{1.,1.},
/*03*/{-1., 1.},/*04*/{ 0.,-1.},/*05*/{1.,0.},
/*06*/{ 0., 1.},/*07*/{-1., 0.},/*08*/{0.,0.} };

int TPZGeoElQ2d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  pss[1] = 0.; pss[2] = 0.;//1d e 2d
  pfs[1] = 0.; pfs[2] = 0.;
  pf[1] = 0.; pf[2] = 0.;
  for(sn=0;sn<4;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<9;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];// = 0
      if(son->WhichSide(ps) != sd) cout << "Lado nao bate\n";
      TPZTransform telsd = TPZShapeQuad::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform sont(son->SideDimension(sd));
      TPZTransform t = son->BuildTransform2(sd,gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = TPZShapeQuad::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(8).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) << endl;
        PZError << "father = " << ((son->Father2(8).Element())->Id()) << endl;
        PZError << "side   = " << sd << endl << endl;
        int ok;
        cin >> ok;
      } else {
      	cout << "Transformacao OK!\n";
       	cout << "Filho/lado : " << son->Id() << "/" << sd << endl;
        cout << "Pai : " << son->Father2(8).Element()->Id() << endl << endl;
      }
    }
  }
  return 1;
}

/*
REAL TPZGeoElQ2d::Mesure(int dim) {
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
    fMesure = fabs(xx + yy);
  }
  return fMesure;
}
*/
//void TPZGeoElQ2d::Center(TPZVec<REAL> &center) {
//  center[0] = center[1] = 0.;
//}
void TPZGeoElQ2d::SetSubElement(int id, TPZGeoEl *el){
  if (id<0 || id >3){
    PZError << "TPZGeoElQ2d::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id]=el;
  return;
}

TPZIntPoints * TPZGeoElQ2d::CreateSideIntegrationRule(int side, int order)
{
	return TPZGeoQuad::CreateSideIntegrationRule(side,order);
}


TPZTransform TPZGeoElQ2d::GetTransform(int side,int son) {
	return TPZRefQuad::GetTransform(side,son);
}

int TPZGeoElQ2d::NSubElements() {
	return TPZRefQuad::NSubEl;
}

void TPZGeoElQ2d::CenterPoint(int side, TPZVec<REAL> &masscent){
  TPZShapeQuad::CenterPoint(side,masscent);
}

REAL TPZGeoElQ2d::RefElVolume(){
  return TPZShapeQuad::RefElVolume();
}
