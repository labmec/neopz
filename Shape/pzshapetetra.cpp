// $Id: pzshapetetra.cpp,v 1.5 2005-02-28 22:11:26 phil Exp $
#include "pzshapetetra.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {

/*Projection of the point within a tetraedro to a rib*/
REAL TPZShapeTetra::gRibTrans3dTetr1d[6][3] = {
  {2.,1.,1.} , {-1.,1.,0.} , {-1.,-2.,-1.} ,
  {1.,1.,2.} , {-1.,0.,1.} , {0.,-1.,1.}
};
REAL TPZShapeTetra::gVet1dTetr[6] = { -1., 0., 1., -1., 0., 0. };
/*Projection of the point within a tetraedro to a face*/
REAL TPZShapeTetra::gFaceTrans3dTetr2d[4][2][3] = {
  { {2.,0.,0.},{0.,2.,0.} },
  { {2.,0.,0.},{0.,0.,2.} },
  { {-2./3.,4./3.,-2./3.},{-2./3.,-2./3.,4./3.} },
  { {0.,2.,0.},{0.,0.,2.} }
};

REAL TPZShapeTetra::gVet2dTetr[4][2] = { {-1.,-1.},{-1.,-1.},{-1./3.,-1./3.},{-1.,-1.} };

REAL TPZShapeTetra::gFaceSum3dTetra2d[4][2] = {
   {-1.,-1.},{-1.,-1.},{-1./3.,-1./3.},{-1.,-1.}
};

REAL TPZShapeTetra::gFaceTrans3dTetra2d[4][2][3] = {
  { { 2.,0.,0.},{0.,2.,0.} },//segundo FaceSides[4][3]
  { { 2.,0.,0.},{0.,0.,2.} },//toma o sentido dos eixos
  { { -2./3.,4./3.,-2./3.},{-2./3.,-2./3.,4./3.} },
  { { 0.,2.,0.},{0.,0.,2.} }
};

REAL TPZShapeTetra::gRibSum3dTetra1d[6] = {-1.,0.,1.,-1.,0.,0.};//estava : {-1.,0.,1.,-1.,0.,0.};

REAL TPZShapeTetra::gRibTrans3dTetra1d[6][3] = {//segundo SideNodes[6][2] vale -1 no 1o extremo e +1 no 2o
  {2.,1.,1.} , {-1.,1.,0. } , {-1.,-2.,-1.} ,//{2.,1.,1.} , {-1.,1.,0. } , {-1.,-2.,-1.} ,
  {1.,1.,2.} , {-1.,0.,1. } , {0.,-1.,1.}//{1.,1.,2.} , {-1.,0.,1. } , {0.,-1.,1.}
};

int TPZShapeTetra::FaceNodes[4][3]  = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

int TPZShapeTetra::SideNodes[6][2]  = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };

int TPZShapeTetra::ShapeFaceId[4][3] = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

static int sidedimension[15] = {0,0,0,0,1,1,1,1,1,1,2,2,2,2,3};

static int nhighdimsides[15] = {7,7,7,7,3,3,3,3,3,3,1,1,1,1,0};

static int highsides[15][7] = {
{4,6,7,10,11,13,14},
{4,5,8,10,11,12,14},
{5,6,9,10,12,13,14},
{7,8,9,11,12,13,14},
{10,11,14},
{10,12,14},
{10,13,14},
{11,13,14},
{11,12,14},
{12,13,14},
{14},
{14},
{14},
{14},
{-999}
};

static REAL sidetosidetransforms[15][7][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
},
{
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,0}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,0}}
},
{
{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,0}}
},
{
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0,0,0.5},{-99,-99,-99},{-99,-99,-99},{0,0,0.5}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-0.5,0,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,0,0.5}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0,0.5,0.5}}
},
{
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,0}}
},
{
{{1,0,0},{0,0,1},{-99,-99,-99},{0,0,0}}
},
{
{{-1,1,0},{-1,0,1},{-99,-99,-99},{1,0,0}}
},
{
{{0,1,0},{0,0,1},{-99,-99,-99},{0,0,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static REAL MidSideNode[15][3] = {
/*00*/{.0,.0},/*01*/{1.,.0},/*02*/{0.,1.,.0},/*03*/{.0,0.,1.0},/*04*/{.5,.0,.0},
/*05*/{.5,.5},/*06*/{0.,.5},/*07*/{0.,0.,.5},/*08*/{.5,0.,0.5},/*09*/{.0,.5,.5},
/*10*/{1./3.,1./3., 0.  }  ,/*11*/{1./3., .0  ,1./3.},
/*12*/{1./3.,1./3.,1./3.}  ,/*13*/{ 0.  ,1./3.,1./3.},/*14*/{1./4.,1./4.,1./4.} };

static int FaceConnectLocId[4][7] = { {0,1,2,4,5,6,10},{0,1,3,4,8,7,11},
					     {1,2,3,5,9,8,12},{0,2,3,6,9,7,13} };

void TPZShapeTetra::CornerShape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  phi(0,0)  = 1-pt[0]-pt[1]-pt[2];
  phi(1,0)  = pt[0];
  phi(2,0)  = pt[1];
  phi(3,0)  = pt[2];

  dphi(0,0) = -1.0;
  dphi(1,0) = -1.0;
  dphi(2,0) = -1.0;
  dphi(0,1) =  1.0;
  dphi(1,1) =  0.0;
  dphi(2,1) =  0.0;
  dphi(0,2) =  0.0;
  dphi(1,2) =  1.0;
  dphi(2,2) =  0.0;
  dphi(0,3) =  0.0;
  dphi(1,3) =  0.0;
  dphi(2,3) =  1.0;
}

//ifstream inn("mats.dt");
void TPZShapeTetra::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {

  CornerShape(pt,phi,dphi);
  //  if(order[9]<2) return;
  int shape = 4;
  //rib shapes
  for (int rib = 0; rib < 6; rib++) {
    REAL outval;
    ProjectPoint3dTetraToRib(rib,pt,outval);
    TPZVec<int> ids(2);
    TPZManVector<REAL,1> outvalvec(1,outval);
    int id0,id1;
    id0 = SideNodes[rib][0];
    id1 = SideNodes[rib][1];
    ids[0] = id[id0];
    ids[1] = id[id1];
    REAL store1[20],store2[60];
    int ordin = order[rib]-1;//three orders : order in x , order in y and order in z
    TPZFMatrix phin(ordin,1,store1,20),dphin(3,ordin,store2,60);
    phin.Zero();
    dphin.Zero();
    TPZShapeLinear::ShapeInternal(outvalvec,order[rib],phin,dphin,TPZShapeLinear::GetTransformId1d(ids));//ordin = ordem de um lado
    TransformDerivativeFromRibToTetra(rib,ordin,dphin);
    for (int i = 0; i < ordin; i++) {
      phi(shape,0) = phi(id0,0)*phi(id1,0)*phin(i,0);
      for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj ,id0) * phi(id1, 0 )  * phin( i, 0) +
                          phi(id0, 0 )  * dphi(xj ,id1) * phin( i, 0) +
                          phi(id0, 0 )  * phi(id1, 0 )  * dphin(xj,i);
      }
      shape++;
    }
  }
  //  if(order[10]<3) return;
  //face shapes
  for (int face = 0; face < 4; face++) {
    if (order[6+face] < 3) continue;
    TPZVec<REAL> outval(2);
    ProjectPoint3dTetraToFace(face,pt,outval);
    REAL store1[20],store2[60];
    int ord1;//,ord2;
    //elt->FaceOrder(face,ord1,ord2);
	ord1 = order[6+face];
	//ord2 = ord1;
    if(ord1<3) continue;
    int ordin =  (ord1-2)*(ord1-1)/2;
    TPZFMatrix phin(ordin,1,store1,20),dphin(3,ordin,store2,60);//ponto na face
    phin.Zero();
    dphin.Zero();
    TPZManVector<int> ids(3);
    int id0,id1,id2,i;
    for(i=0;i<3;i++) ids[i] = id[FaceNodes[face][i]];
    id0 = ShapeFaceId[face][0];//indice das shapes da face que compoem a shape atual
    id1 = ShapeFaceId[face][1];//equivale a FaceIdsCube(face,ids,id,id0,id1);
    id2 = ShapeFaceId[face][2];
    int transid = TPZShapeTriang::GetTransformId2dT(ids);
    outval[0] = (outval[0]+1.)/2.;//devido a correção na função
    outval[1] = (outval[1]+1.)/2.;//Shape2dTriangleInternal(..)
    TPZShapeTriang::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
    int c = dphin.Cols();
    for(i=0;i<c;i++) {
    	dphin(0,i) /= 2.;
      dphin(1,i) /= 2.;
      dphin(2,i) /= 2.;
    }
    TransformDerivativeFromFaceToTetra(face,ordin,dphin);//ord = numero de shapes
    for(i=0;i<ordin;i++)	{
      phi(shape,0) = phi(id0,0)*phi(id1,0)*phi(id2,0)*phin(i,0);
      for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj,id0)* phi(id1 , 0 )* phi(id2 , 0 )* phin(i ,0) +
                           phi(id0, 0)*dphi(xj  ,id1)* phi(id2 , 0 )* phin(i ,0) +
                           phi(id0, 0)* phi(id1 , 0 )*dphi(xj  ,id2)* phin(i ,0) +
                           phi(id0, 0)* phi(id1 , 0 )* phi(id2 , 0 )*dphin(xj,i);
      }
      shape++;
    }
  }
  if(order[10]<4) return;
  //volume shapes
  REAL store1[20],store2[60];
  int totsum = 0,sum;
  int i;
  for(i=0;i<order[10]-3;i++) {
    sum = (i+1)*(i+2) / 2;
    totsum += sum;
  }
  int ord = totsum;
  TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);
  phin.Zero();
  dphin.Zero();
  ShapeInternal(pt,order[10],phin,dphin);
  for(i=0;i<ord;i++)	{
    phi(shape,0) = phi(0,0)*phi(1,0)*phi(2,0)*phi(3,0)*phin(i,0);
    for(int xj=0;xj<3;xj++) {
      dphi(xj,shape) = dphi(xj,0)* phi(1 ,0)* phi(2 ,0)* phi(3 ,0)* phin(i ,0) +
                        phi(0, 0)*dphi(xj,1)* phi(2 ,0)* phi(3 ,0)* phin(i ,0) +
                        phi(0, 0)* phi(1 ,0)*dphi(xj,2)* phi(3 ,0)* phin(i ,0) +
                        phi(0, 0)* phi(1 ,0)* phi(2 ,0)*dphi(xj,3)* phin(i ,0) +
                        phi(0, 0)* phi(1 ,0)* phi(2 ,0)* phi(3 ,0)*dphin(xj,i);
    }
    shape++;
  }
}

void TPZShapeTetra::ProjectPoint3dTetrSide(int side, TPZVec<REAL> &in, REAL &out) {

  out = gRibTrans3dTetr1d[side][0]*in[0]+gRibTrans3dTetr1d[side][1]*in[1]+gRibTrans3dTetr1d[side][2]*in[2]+gVet1dTetr[side];
}

void TPZShapeTetra::ProjectPoint3dTetrFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out) {

  out[0] = gFaceTrans3dTetr2d[face][0][0]*in[0]+gFaceTrans3dTetr2d[face][0][1]*in[1]+gFaceTrans3dTetr2d[face][0][2]*in[2]+gVet2dTetr[face][0];
  out[1] = gFaceTrans3dTetr2d[face][1][0]*in[0]+gFaceTrans3dTetr2d[face][1][1]*in[1]+gFaceTrans3dTetr2d[face][1][2]*in[2]+gVet2dTetr[face][1];
}

void TPZShapeTetra::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                    TPZFMatrix &dphi) {
  if(order < 4) return;
  int ord = order-3;
  REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
  TPZFMatrix phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
            dphi0(1,ord,store4,20),dphi1(1,ord,store5,20),dphi2(1,ord,store6,20);
  TPZShapeLinear::fOrthogonal(2.*x[0]-1.,ord,phi0,dphi0);
  TPZShapeLinear::fOrthogonal(2.*x[1]-1.,ord,phi1,dphi1);
  TPZShapeLinear::fOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);
  int index = 0;
  for (int i=0;i<ord;i++) {
    for (int j=0;j<ord;j++) {
      for (int k=0;k<ord;k++) {
      	if( i+j+k < ord ) {
            //int index = ord*(ord*i+j)+k;
             phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
            dphi(0,index) =  2.*dphi0(0,i)* phi1(j,0)* phi2(k,0);
            dphi(1,index) =  2.* phi0(i,0)*dphi1(0,j)* phi2(k,0);
            dphi(2,index) =  2.* phi0(i,0)* phi1(j,0)*dphi2(0,k);
            index++;
         }
      }
    }
  }
}

void TPZShapeTetra::TransformDerivativeFromRibToTetra(int rib,int num,TPZFMatrix &dphi) {
  for (int j = 0;j<num;j++) {
    dphi(2,j) = gRibTrans3dTetra1d[rib][2]*dphi(0,j);
    dphi(1,j) = gRibTrans3dTetra1d[rib][1]*dphi(0,j);
    dphi(0,j) = gRibTrans3dTetra1d[rib][0]*dphi(0,j);
  }
}

void TPZShapeTetra::TransformDerivativeFromFaceToTetra(int face,int num,TPZFMatrix &dphi) {

  for (int j = 0;j<num;j++) {
    dphi(2,j) = gFaceTrans3dTetra2d[face][0][2]*dphi(0,j)+gFaceTrans3dTetra2d[face][1][2]*dphi(1,j);
    REAL dphi1j = dphi(1,j);
    dphi(1,j) = gFaceTrans3dTetra2d[face][0][1]*dphi(0,j)+gFaceTrans3dTetra2d[face][1][1]*dphi(1,j);
    dphi(0,j) = gFaceTrans3dTetra2d[face][0][0]*dphi(0,j)+gFaceTrans3dTetra2d[face][1][0]*dphi1j;//dphi(1,j);
  }
}

void TPZShapeTetra::ProjectPoint3dTetraToRib(int rib, TPZVec<REAL> &in, REAL &outval) {
  outval = gRibTrans3dTetra1d[rib][0]*in[0]+gRibTrans3dTetra1d[rib][1]*in[1]+gRibTrans3dTetra1d[rib][2]*in[2]+gRibSum3dTetra1d[rib];
}

void TPZShapeTetra::ProjectPoint3dTetraToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval) {
  outval[0] = gFaceTrans3dTetra2d[face][0][0]*in[0]+gFaceTrans3dTetra2d[face][0][1]*in[1]+gFaceTrans3dTetra2d[face][0][2]*in[2]+gFaceSum3dTetra2d[face][0];
  outval[1] = gFaceTrans3dTetra2d[face][1][0]*in[0]+gFaceTrans3dTetra2d[face][1][1]*in[1]+gFaceTrans3dTetra2d[face][1][2]*in[2]+gFaceSum3dTetra2d[face][1];
}

int TPZShapeTetra::NConnectShapeF(int side, int order) {
   if(side<4) return 1;//0 a 3
   //   int s = side-4;
   if(side<10) return order-1;//4 a 9
   if(side<14) {//10 a 13
   	int sum = 0;
      for(int i=0;i<order-1;i++) sum += i;
   	return sum;
   }
   if(side==14) {
   	int totsum = 0,sum;
      for(int i=1;i<order-2;i++) {
         sum = i*(i+1) / 2;
         totsum += sum;
      }
      return totsum;
   }
   PZError << "TPZCompElT3d::NConnectShapeF, bad parameter side " << side << endl;
   return 0;
}

int TPZShapeTetra::NShapeF(TPZVec<int> &order) {
  int in,res=NNodes;
  for(in=NNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NNodes]);
  return res;
}

int TPZShapeTetra::NConnects() {
	return NSides;
}

TPZTransform TPZShapeTetra::TransformElementToSide(int side){

	if(side<0 || side>14){
  	PZError << "TPZShapeTetra::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

  TPZTransform t(sidedimension[side],3);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
    case 0:
    case 1:
    case 2:
    case 3:
      return t;
    case 4:
    	t.Mult()(0,0) =  2.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 5:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) =  1.0;
      return t;
    case 6:
    	t.Mult()(0,1) = -2.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case 7:
    	t.Mult()(0,2) =  2.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 8:
    	t.Mult()(0,0) = -1.0;
      t.Mult()(0,2) =  1.0;
      return t;
    case 9:
      t.Mult()(0,1) = -1.0;
      t.Mult()(0,2) =  1.0;
      return t;
    case 10:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
    case 11:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,2) =  1.0;
      return t;
    case 12:
    case 13:
      t.Mult()(0,1) =  1.0;
      t.Mult()(1,2) =  1.0;
      return t;
    case 14:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapeTetra::TransformSideToElement(int side){

	if(side<0 || side>14){
  	PZError << "TPZShapeTetra::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(3,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      return t;
    case 2:
      t.Sum()(1,0) =  1.0;
      return t;
    case 3:
    	t.Sum()(2,0) =  1.0;
      return t;
    case 4:
      t.Mult()(0,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      return t;
    case 5:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 6:
      t.Mult()(1,0) = -0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 7:
      t.Mult()(2,0) =  0.5;
      t.Sum() (2,0) =  0.5;
      return t;
    case 8:
      t.Mult()(0,0) = -0.5;
      t.Mult()(2,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (2,0) =  0.5;
      return t;
    case 9:
      t.Mult()(1,0) = -0.5;
      t.Mult()(2,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      t.Sum() (2,0) =  0.5;
      return t;
    case 10:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
    case 11:
      t.Mult()(0,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      return t;
    case 12:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) = -1.0;
      t.Mult()(1,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      t.Sum() (0,0) =  1.0;
      return t;
    case 13:
      t.Mult()(1,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      return t;
    case 14:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;

  }
	return TPZTransform(0,0);
}


int TPZShapeTetra::SideNodeLocId(int side, int node)
{
	if(side<4 && node == 0) return side;
	if(side>=4 && side < 10 && node <2) return SideNodes[side-4][node];
	if(side >= 10 && side < 14 && node <3) return FaceNodes[side-10][node];
	if(side ==14 && node < 4) return node;
	PZError << "TPZShapeTetra::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;

}

static int nsidenodes[15] = 
{
	1,1,1,1,
	2,2,2,2,2,2,
	3,3,3,3,
	4};
int TPZShapeTetra::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapeTetra::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapeTetra::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapeTetra::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapeTetra::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
		return TPZTransform(0);
	}
	if(sidefrom == sideto) {
		return TPZTransform(sidedimension[sidefrom]);
	}
	if(sidefrom == NSides-1) {
		return TransformElementToSide(sideto);
	}
	int nhigh = nhighdimsides[sidefrom];
	int is;
	for(is=0; is<nhigh; is++) {
		if(highsides[sidefrom][is] == sideto) {
			int dfr = sidedimension[sidefrom];
			int dto = sidedimension[sideto];
			TPZTransform trans(dto,dfr);
			int i,j;
			for(i=0; i<dto; i++) {
				for(j=0; j<dfr; j++) {
					trans.Mult()(i,j) = sidetosidetransforms[sidefrom][is][j][i];
				}
				trans.Sum()(i,0) = sidetosidetransforms[sidefrom][is][3][i];
			}
			return trans;
		}
	}
	PZError << "TPZShapeTetra::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapeTetra::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapeTetra::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

int TPZShapeTetra::NSideConnects(int side) {
	if(side<0)   return -1;
	if(side<4)   return 1;//cantos : 0 a 3
   if(side<10)  return 3;//lados : 4 a 9
   if(side<14)	 return 7;//faces : 10 a 13
   if(side==14) return 15;//centro : 15
   return -1;
}

int TPZShapeTetra::SideConnectLocId(int side, int node) {
	if(side<0 || side>15) return -1;
   if(side<4) {
   	if(node==0) return side;
   } else
   if(side<7) {//4,5,6
	   int s = side-4;//0,1,2
   	if(!node) return s;//0,1,2
      if(node==1) return (s+1)%3;//1,2,0
      if(node==2) return side;//4,5,6
   } else
   if(side<10) {//7,8,9
   	int s = side-7;//0,1,2
   	if(!node) return s;//0,1,2
      if(node==1) return 3;//4,4,4
      if(node==2) return side;//7,8,9
   } else
   if(side<14) {//10 a 13
   	int s = side-10;
   	if(node<7) return FaceConnectLocId[s][node];
   } else
   if(side==14 && node<15){
   	return node;
   }
	PZError << "TPZShapeTetra::SideConnectLocId called for node = "
   	      << node << " and side = " << side << "\n";
   return -1;
}

void TPZShapeTetra::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}

void TPZShapeTetra::SideShape(int side, TPZVec<REAL> &point, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {

  if(side<0 || side>15) PZError << "TPZCompElT3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==14) Shape(point,id,order,phi,dphi);
  else if(side<4) phi(0,0)=1.;
  else if(side<10) {//4 a 9
    TPZShapeLinear::Shape(point,id,order,phi,dphi);
  }
  else if(side<14) {//faces 10,11,12,13
    TPZShapeTriang::Shape(point,id,order,phi,dphi);
  }

}

};
