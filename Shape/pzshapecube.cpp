#include "pzshapecube.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

namespace pzshape {

REAL TPZShapeCube::gRibTrans3dCube1d[12][3] = {
  {1., 0.,0.} , { 0.,1.,0.} , {-1., 0.,0.} ,
  {0.,-1.,0.} , { 0.,0.,1.} , { 0., 0.,1.} ,
  {0., 0.,1.} , { 0.,0.,1.} , { 1., 0.,0.} ,
  {0., 1.,0.} , {-1.,0.,0.} , { 0.,-1.,0.}
};

REAL TPZShapeCube::gFaceTrans3dCube2d[6][2][3] = {
  { { 1.,0.,0.},{0.,1.,0.} },
  { { 1.,0.,0.},{0.,0.,1.} },
  { { 0.,1.,0.},{0.,0.,1.} },
  { { 1.,0.,0.},{0.,0.,1.} },//{-1.,0.,0.},{0.,0.,1.}
  { { 0.,1.,0.},{0.,0.,1.} },
  { { 1.,0.,0.},{0.,1.,0.} }
};

int TPZShapeCube::FaceNodes[6][4]  = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7} };

int TPZShapeCube::SideNodes[12][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},
                      	               {2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };


int TPZShapeCube::ShapeFaceId[6][2] = { {0,2},{0,5},{1,6},{3,6},{0,7},{4,6} };

static int sidedimension[27] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3};

static int nhighdimsides[27] = {7,7,7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,0};

static int highsides[27][7] = {
{8,11,12,20,21,24,26},
{8,9,13,20,21,22,26},
{9,10,14,20,22,23,26},
{10,11,15,20,23,24,26},
{12,16,19,21,24,25,26},
{13,16,17,21,22,25,26},
{14,17,18,22,23,25,26},
{15,18,19,23,24,25,26},
{20,21,26},
{20,22,26},
{20,23,26},
{20,24,26},
{21,24,26},
{21,22,26},
{22,23,26},
{23,24,26},
{21,25,26},
{22,25,26},
{23,25,26},
{24,25,26},
{26},
{26},
{26},
{26},
{26},
{26},
{-999}
};

static REAL sidetosidetransforms[27][7][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,1}}
},
{
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,-1}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,-1}}
},
{
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,-1}}
},
{
{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,-1}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{-1,-1,0}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,-1,0}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,1,0}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{-1,1,0}}
},
{
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,1}}
},
{
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,1}}
},
{
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,1}}
},
{
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,1}}
},
{
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,-1}}
},
{
{{1,0,0},{0,0,1},{-99,-99,-99},{0,-1,0}}
},
{
{{0,1,0},{0,0,1},{-99,-99,-99},{1,0,0}}
},
{
{{1,0,0},{0,0,1},{-99,-99,-99},{0,1,0}}
},
{
{{0,1,0},{0,0,1},{-99,-99,-99},{-1,0,0}}
},
{
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,1}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static int FaceSides[6][4] = { {8,9,10,11},{8,13,16,12},{9,14,17,13},
				      {10,14,18,15},{11,15,19,12},{16,17,18,19} };

static int FaceConnectLocId[6][9] = { {0,1,2,3,8,9,10,11,20},{0,1,5,4,8,13,16,12,21},
					     {1,2,6,5,9,14,17,13,22},{3,2,6,7,10,14,18,15,23},//{2,3,7,6,10,15,18,14,23}
					     {0,3,7,4,11,15,19,12,24},{4,5,6,7,16,17,18,19,25} };

static REAL MidSideNode[27][3] = {
/*00*/{-1.,-1.,-1.},/*01*/{1.,-1.,-1.},/*02*/{1.,1.,-1.},/*03*/{-1.,1.,-1.},
/*04*/{-1.,-1., 1.},/*05*/{1.,-1., 1.},/*06*/{1.,1., 1.},/*07*/{-1.,1., 1.},
/*08*/{ 0.,-1.,-1.},/*09*/{1., 0.,-1.},/*10*/{0.,1.,-1.},/*11*/{-1.,0.,-1.},
/*12*/{-1.,-1., 0.},/*13*/{1.,-1., 0.},/*14*/{1.,1., 0.},/*15*/{-1.,1., 0.},
/*16*/{ 0.,-1., 1.},/*17*/{1., 0., 1.},/*18*/{0.,1., 1.},/*19*/{-1.,0., 1.},
/*20*/{ 0., 0.,-1.},/*21*/{0.,-1., 0.},/*22*/{1.,0., 0.},/*23*/{ 0.,1., 0.},
/*24*/{-1., 0., 0.},/*25*/{0., 0., 1.},/*26*/{0.,0., 0.} };

void TPZShapeCube::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {

      REAL x[2],dx[2],y[2],dy[2],z[2],dz[2];
      x[0]  = (1.-pt[0])/2.;
      x[1]  = (1.+pt[0])/2.;
      dx[0] = -0.5;
      dx[1] =  0.5;
      y[0]  = (1.-pt[1])/2.;
      y[1]  = (1.+pt[1])/2.;
      dy[0] = -0.5;
      dy[1] =  0.5;
      z[0]  = (1.-pt[2])/2.;
      z[1]  = (1.+pt[2])/2.;
      dz[0] = -0.5;
      dz[1] =  0.5;

      phi(0,0)  = x[0]*y[0]*z[0];
      phi(1,0)  = x[1]*y[0]*z[0];
      phi(2,0)  = x[1]*y[1]*z[0];
      phi(3,0)  = x[0]*y[1]*z[0];
      phi(4,0)  = x[0]*y[0]*z[1];
      phi(5,0)  = x[1]*y[0]*z[1];
      phi(6,0)  = x[1]*y[1]*z[1];
      phi(7,0)  = x[0]*y[1]*z[1];
      dphi(0,0) = dx[0]*y[0]*z[0];
      dphi(1,0) = x[0]*dy[0]*z[0];
      dphi(2,0) = x[0]*y[0]*dz[0];
      dphi(0,1) = dx[1]*y[0]*z[0];
      dphi(1,1) = x[1]*dy[0]*z[0];
      dphi(2,1) = x[1]*y[0]*dz[0];
      dphi(0,2) = dx[1]*y[1]*z[0];
      dphi(1,2) = x[1]*dy[1]*z[0];
      dphi(2,2) = x[1]*y[1]*dz[0];
      dphi(0,3) = dx[0]*y[1]*z[0];
      dphi(1,3) = x[0]*dy[1]*z[0];
      dphi(2,3) = x[0]*y[1]*dz[0];
      dphi(0,4) = dx[0]*y[0]*z[1];
      dphi(1,4) = x[0]*dy[0]*z[1];
      dphi(2,4) = x[0]*y[0]*dz[1];
      dphi(0,5) = dx[1]*y[0]*z[1];
      dphi(1,5) = x[1]*dy[0]*z[1];
      dphi(2,5) = x[1]*y[0]*dz[1];
      dphi(0,6) = dx[1]*y[1]*z[1];
      dphi(1,6) = x[1]*dy[1]*z[1];
      dphi(2,6) = x[1]*y[1]*dz[1];
      dphi(0,7) = dx[0]*y[1]*z[1];
      dphi(1,7) = x[0]*dy[1]*z[1];
      dphi(2,7) = x[0]*y[1]*dz[1];
}


void TPZShapeCube::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {
  ShapeCorner(pt,phi,dphi);
  int shape = 8;
  //rib shapes
  for (int rib = 0; rib < 12; rib++) {
    REAL outval;
    ProjectPoint3dCubeToRib(rib,pt,outval);
    TPZManVector<REAL,1> outvalvec(1,outval);
    TPZVec<int> ids(2);
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
    TransformDerivativeFromRibToCube(rib,ordin,dphin);
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
  //face shapes
  for (int face = 0; face < 6; face++) {

    TPZVec<REAL> outval(2);
    ProjectPoint3dCubeToFace(face,pt,outval);
    REAL store1[20],store2[60];
    int ord1,ord2;
	ord1 = order[12+face];
	ord2=ord1;
//    FaceOrder(face,ord1,ord2);
    if(ord1<2 || ord2<2) continue;
    int ord =  (ord1-1)*(ord2-1);
    TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);//ponto na face
    phin.Zero();
    dphin.Zero();
    int ordin =  (ord1 > ord2) ? ord1 : ord2;
    ordin--;
    TPZManVector<int> ids(4);
    //TPZVec<int> ids(4);
    int id0,id1,i;
    for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
    id0 = ShapeFaceId[face][0];//numero das shapes da face que compoem a shape atual
    id1 = ShapeFaceId[face][1];
    TPZShapeQuad::ShapeInternal(outval,ord1-2,phin,dphin,TPZShapeQuad::GetTransformId2dQ(ids));//ordin = ordem de um lado
    TransformDerivativeFromFaceToCube(face,ord,dphin);//ord = numero de shapes
    for(i=0;i<ord;i++)	{
      phi(shape,0) = phi(id0,0)*phi(id1,0)*phin(i,0);
      for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj,id0)* phi(id1 , 0 )* phin(i ,0) +
                           phi(id0, 0)*dphi(xj  ,id1)* phin(i ,0) +
                           phi(id0, 0)* phi(id1 , 0 )*dphin(xj,i);
      }
      shape++;
    }
  }
  //volume shapes
  REAL store1[20],store2[60];
  int ordmin1 = (order[18]-1);
  int ord =  ordmin1*ordmin1*ordmin1;//(p-1)^3 : 0<=n1,n2,n3<=p-2
  TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);
  phin.Zero();
  dphin.Zero();
  ShapeInternal(pt,ordmin1,phin,dphin);
  for(int i=0;i<ord;i++)	{
    phi(shape,0) = phi(0,0)*phi(6,0)*phin(i,0);
    for(int xj=0;xj<3;xj++) {
      dphi(xj,shape) = dphi(xj,0)* phi(6 ,0)* phin(i ,0) +
                        phi(0, 0)*dphi(xj,6)* phin(i ,0) +
                        phi(0, 0)* phi(6 ,0)*dphin(xj,i);
    }
    shape++;
  }
}
void TPZShapeCube::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {
  if(side<0 || side>26) PZError << "TPZCompElC3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==26) Shape(pt,id,order,phi,dphi);
  else if(side<8) TPZShapePoint::Shape(pt,id,order,phi,dphi);
  else if(side<20) {//8 a 19
    TPZShapeLinear::Shape(pt,id,order,phi,dphi);
  }
  else if(side<26) {//faces
    TPZShapeQuad::Shape(pt,id,order,phi,dphi);
  }
}

void TPZShapeCube::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                    TPZFMatrix &dphi) {//,int cube_transformation_index
  if(order < 1) return;
  int ord = order;//fSideOrder[18]-1;
  order = order*order*order;
  phi.Resize(order,1);
  dphi.Resize(3,order);
  REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
  TPZFMatrix phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
    dphi0(1,ord,store4,20),dphi1(1,ord,store5,20),dphi2(1,ord,store6,20);
  TPZShapeLinear::fOrthogonal(x[0],ord,phi0,dphi0);
  TPZShapeLinear::fOrthogonal(x[1],ord,phi1,dphi1);
  TPZShapeLinear::fOrthogonal(x[2],ord,phi2,dphi2);
  for (int i=0;i<ord;i++) {
    for (int j=0;j<ord;j++) {
      for (int k=0;k<ord;k++) {
         int index = ord*(ord*i+j)+k;
          phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
         dphi(0,index) = dphi0(0,i)* phi1(j,0)* phi2(k,0);
         dphi(1,index) =  phi0(i,0)*dphi1(0,j)* phi2(k,0);
         dphi(2,index) =  phi0(i,0)* phi1(j,0)*dphi2(0,k);
      }
    }
  }
}

void TPZShapeCube::TransformDerivativeFromRibToCube(int rib,int num,TPZFMatrix &dphi) {
  for (int j = 0;j<num;j++) {
    dphi(2,j) = gRibTrans3dCube1d[rib][2]*dphi(0,j);
    dphi(1,j) = gRibTrans3dCube1d[rib][1]*dphi(0,j);
    dphi(0,j) = gRibTrans3dCube1d[rib][0]*dphi(0,j);
  }
}

void TPZShapeCube::ProjectPoint3dCubeToRib(int side, TPZVec<REAL> &in, REAL &outval) {
  outval = gRibTrans3dCube1d[side][0]*in[0]+gRibTrans3dCube1d[side][1]*in[1]+gRibTrans3dCube1d[side][2]*in[2];
}

void TPZShapeCube::ProjectPoint3dCubeToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval) {
  outval[0] = gFaceTrans3dCube2d[face][0][0]*in[0]+gFaceTrans3dCube2d[face][0][1]*in[1]+gFaceTrans3dCube2d[face][0][2]*in[2];
  outval[1] = gFaceTrans3dCube2d[face][1][0]*in[0]+gFaceTrans3dCube2d[face][1][1]*in[1]+gFaceTrans3dCube2d[face][1][2]*in[2];
}

void TPZShapeCube::TransformDerivativeFromFaceToCube(int rib,int num,TPZFMatrix &dphi) {

  for (int j = 0;j<num;j++) {
    dphi(2,j) = gFaceTrans3dCube2d[rib][0][2]*dphi(0,j)+gFaceTrans3dCube2d[rib][1][2]*dphi(1,j);
    REAL dphi1j = dphi(1,j);
    dphi(1,j) = gFaceTrans3dCube2d[rib][0][1]*dphi(0,j)+gFaceTrans3dCube2d[rib][1][1]*dphi(1,j);
    dphi(0,j) = gFaceTrans3dCube2d[rib][0][0]*dphi(0,j)+gFaceTrans3dCube2d[rib][1][0]*dphi1j;//dphi(1,j);
  }
}

void TPZShapeCube::ProjectPoint3dCubeSide(int side, TPZVec<REAL> &in, REAL &out) {

  out = gRibTrans3dCube1d[side][0]*in[0]+gRibTrans3dCube1d[side][1]*in[1]+gRibTrans3dCube1d[side][2]*in[2];
}

void TPZShapeCube::ProjectPoint3dCubeFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out) {

  out[0] = gFaceTrans3dCube2d[face][0][0]*in[0]+gFaceTrans3dCube2d[face][0][1]*in[1]+gFaceTrans3dCube2d[face][0][2]*in[2];
  out[1] = gFaceTrans3dCube2d[face][1][0]*in[0]+gFaceTrans3dCube2d[face][1][1]*in[1]+gFaceTrans3dCube2d[face][1][2]*in[2];
}

int TPZShapeCube::NConnectShapeF(int side, int order){
   if(side<8) return 1;//0 a 4
   if(side<20) return (order-1);//6 a 14
   if(side<26) {
      return ((order-1)*(order-1));
   }
   if(side==26) {
      return ((order-1)*(order-1)*(order-1));
   }
   PZError << "TPZShapeCube::NConnectShapeF, bad parameter side " << side << endl;
   return 0;
}

int TPZShapeCube::NShapeF(TPZVec<int> &order) {
  int in,res=NNodes;
  for(in=NNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NNodes]);
  return res;
}

int TPZShapeCube::NConnects() {
	return 27;
}

TPZTransform TPZShapeCube::TransformElementToSide(int side){

	if(side<0 || side>26){
  	PZError << "TPZShapeCube::TransformElementToSide called with side error\n";
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
    case 4:
    case 5 :
    case 6:
    case 7:
      return t;
    case  8:
    case 16:
    	t.Mult()(0,0) =  1.0;
      return t;
    case  9:
    case 17:
    	t.Mult()(0,1) =  1.0;
      return t;
    case 10:
    case 18:
    	t.Mult()(0,0) = -1.0;
      return t;
    case 11:
    case 19:
    	t.Mult()(0,1) = -1.0;
      return t;
    case 12:
    case 13:
    case 14:
    case 15:
    	t.Mult()(0,2) = 1.0;
      return t;
    case 20:
    case 25:
    	t.Mult()(0,0) =  1.0;
    	t.Mult()(1,1) =  1.0;
      return t;
    case 21:
    case 23:
    	t.Mult()(0,0) =  1.0;
    	t.Mult()(1,2) =  1.0;
      return t;
    case 22:
    case 24:
    	t.Mult()(0,1) =  1.0;
    	t.Mult()(1,2) =  1.0;
      return t;
    case 26:
    	t.Mult()(0,0) =  1.0;
    	t.Mult()(1,1) =  1.0;
    	t.Mult()(2,2) =  1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapeCube::TransformSideToElement(int side){

	if(side<0 || side>26){
  	PZError << "TPZShapeCube::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(3,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) = -1.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) = -1.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 2:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 3:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 4:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) = -1.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 5:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) = -1.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 6:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 7:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 8:
      t.Mult()(0,0) =  1.0;
      t.Sum()(1,0)  = -1.0;
      t.Sum()(2,0)  = -1.0;
      return t;
    case 9:
      t.Mult()(1,0) =  1.0;
      t.Sum()(0,0)  =  1.0;
      t.Sum()(2,0)  = -1.0;
      return t;
    case 10:
      t.Mult()(0,0) = -1.0;
      t.Sum()(1,0)  =  1.0;
      t.Sum()(2,0)  = -1.0;
      return t;
    case 11:
      t.Mult()(1,0) = -1.0;
      t.Sum()(0,0)  = -1.0;
      t.Sum()(2,0)  = -1.0;
      return t;
    case 12:
      t.Mult()(2,0) =  1.0;
      t.Sum()(0,0)  = -1.0;
      t.Sum()(1,0)  = -1.0;
      return t;
    case 13:
      t.Mult()(2,0) =  1.0;
      t.Sum()(0,0)  =  1.0;
      t.Sum()(1,0)  = -1.0;
      return t;
    case 14:
      t.Mult()(2,0) =  1.0;
      t.Sum()(0,0)  =  1.0;
      t.Sum()(1,0)  =  1.0;
      return t;
    case 15:
      t.Mult()(2,0) =  1.0;
      t.Sum()(0,0)  = -1.0;
      t.Sum()(1,0)  =  1.0;
      return t;
    case 16:
      t.Mult()(0,0) =  1.0;
      t.Sum()(1,0)  = -1.0;
      t.Sum()(2,0)  =  1.0;
      return t;
    case 17:
      t.Mult()(1,0) =  1.0;
      t.Sum()(0,0)  =  1.0;
      t.Sum()(2,0)  =  1.0;
      return t;
    case 18:
      t.Mult()(0,0) = -1.0;
      t.Sum()(1,0)  =  1.0;
      t.Sum()(2,0)  =  1.0;
      return t;
    case 19:
      t.Mult()(1,0) = -1.0;
      t.Sum()(0,0)  = -1.0;
      t.Sum()(2,0)  =  1.0;
      return t;
    case 20:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Sum()(2,0)  = -1.0;
      return t;
    case 21:
      t.Mult()(0,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      t.Sum()(1,0)  = -1.0;
      return t;
    case 22:
      t.Mult()(1,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case 23:
      t.Mult()(0,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      t.Sum()(1,0)  =  1.0;
      return t;
    case 24:
      t.Mult()(1,0) =  1.0;
      t.Mult()(2,1) =  1.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 25:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Sum()(2,0)  =  1.0;
      return t;
    case 26:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}

int TPZShapeCube::SideNodeLocId(int side, int node)
{
	if(side<8 && node == 0) return side;
	if(side>=8 && side < 20 && node < 2) return SideNodes[side-8][node];
	if(side>=20 && side < 26 && node < 4) return FaceNodes[side-20][node];
	if(side == 26 && node < 8) return node;
	PZError << "TPZShapeCube::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

static int nsidenodes[27] = {1,1,1,1,1,1,1,1,
2,2,2,2,2,2,2,2,2,2,2,2,
4,4,4,4,4,4,
8};

int TPZShapeCube::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapeCube::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapeCube::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapeCube::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapeCube::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZShapeCube::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapeCube::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapeCube::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}



void TPZShapeCube::LowerDimensionSides(int side,TPZStack<int> &smallsides) {

	cout << "TPZShapeCube::LowerDimensionSides Nao deve ser usado";
	exit(-1);

  if (side < 8) return;
  int i;
  if(side < 20) {//side = 8 a 19 : entram os cantos dos lados
    int s = side-8;
    smallsides.Push(SideNodes[s][0]);
    smallsides.Push(SideNodes[s][1]);
  } else if(side < 26) {//side = 20 a 25
    int s = side-20;   //entram cantos e lados da face
    smallsides.Push(FaceNodes[s][0]);
    smallsides.Push(FaceNodes[s][1]);
    smallsides.Push(FaceNodes[s][2]);
    smallsides.Push(FaceNodes[s][3]);
    smallsides.Push(FaceSides[s][0]);
    smallsides.Push(FaceSides[s][1]);
    smallsides.Push(FaceSides[s][2]);
    smallsides.Push(FaceSides[s][3]);
    //    smallsides.Push(TPZGeoElSide(geo,side));
  } else if(side==26) {//entram todos os cantos, arestas e faces
    for (i=0;i<25;i++) smallsides.Push(i);
  }
}

int TPZShapeCube::NSideConnects(int side) {
  if(side<0) return -1;
  if(side<8) return 1;//cantos : 0 a 7
  if(side<20)	return 3;//lados : 8 a 19
  if(side<26)	return 9;//faces : 20 a 25
  if(side==26)	return 27;//centro : 26
  return -1;
}

// Pronto 23/04/98
int TPZShapeCube::SideConnectLocId(int side, int node) {
  if(side<0 || side>26) return -1;
  if(side<8) {
    if(node==0) return side;
  } 
  else if(side<12) {//8,9,10,11
    int s = side-8;//0,1,2,3
    if(!node) return s;
    if(node==1) return (s+1)%4;
    if(node==2) return side;
  } 
  else if(side<16) {//12,13,14,15
    int s = side-12;//0,1,2,3
    if(!node) return s;
    if(node==1) return s+4;
    if(node==2) return side;
  } 
  else if(side<20) {//16,17,18,19
    int s = side-16;//0,1,2,3
    if(!node) return s+4;
    if(node==1) return (s+1)%4+4;
    if(node==2) return side;
  } 
  else if(side<26) {//20 a 25
    int s = side-20;
    if(node<9) return FaceConnectLocId[s][node];
  } 
  else if(side==26){
    return node;
  }
  PZError << "TPZShapeCube::SideConnectLocId called for node = "
	  << node << " and side = " << side << "\n";
  return -1;
}

void TPZShapeCube::CenterPoint(int side, TPZVec<REAL> &center) {
  center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}


#ifdef _AUTODIFF

void TPZShapeCube::ShapeCube(TPZVec<REAL> &point, TPZVec<int> &id, TPZVec<int> &order, TPZVec<FADREAL> &phi)
{
  const int ndim = 3;

  TPZVec<FADREAL> pt(3);
  pt[0] = point[0];
  pt[0].diff(0, ndim);

  pt[1] = point[1];
  pt[1].diff(1, ndim);

  pt[2] = point[2];
  pt[2].diff(2, ndim);

  ShapeCornerCube(pt,phi);

  int shape = 8;
  //rib shapes
  for (int rib = 0; rib < 12; rib++) {
    FADREAL outval(ndim, 0.0);
    ProjectPoint3dCubeToRib(rib,pt,outval);
    TPZVec<int> ids(2);
    int id0,id1;
    id0 = SideNodes[rib][0];
    id1 = SideNodes[rib][1];
    ids[0] = id[id0];
    ids[1] = id[id1];
    //REAL store1[20], store2[60];
    int ordin = order[rib]-1;//three orders : order in x , order in y and order in z
    //TPZFMatrix phin(ordin,1,store1,20),dphin(3,ordin,store2,60);
    //phin.Zero();
    //dphin.Zero();
    TPZVec<FADREAL> phin(20, FADREAL(ndim, 0.0)); //3d
    TPZShapeLinear::ShapeInternal(outval,ordin,phin,TPZShapeLinear::GetTransformId1d(ids));//ordin = ordem de um lado
//    TransformDerivativeFromRibToCube(rib,ordin,phin);
    for (int i = 0; i < ordin; i++) {
      //phi(shape,0) = phi(id0,0)*phi(id1,0)*phin(i,0);
      phi[shape] = phi[id0] * phi[id1] * phin[i];
      /*for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj ,id0) * phi(id1, 0 )  * phin( i, 0) +
                          phi(id0, 0 )  * dphi(xj ,id1) * phin( i, 0) +
                          phi(id0, 0 )  * phi(id1, 0 )  * dphin(xj,i);
      }*/ // implicitly done
      shape++;
    }

  }
  //face shapes
  for (int face = 0; face < 6; face++) {

    //TPZVec<REAL> outval(2);
    TPZVec<FADREAL> outval(2, FADREAL(ndim, 0.0));
    ProjectPoint3dCubeToFace(face,pt,outval);
  //  REAL store1[20],store2[60];
    int ord1,ord2;
	ord1 = order[12+face];
	ord2=ord1;
//    FaceOrder(face,ord1,ord2); // already commented in the non-FAD version
    if(ord1<2 || ord2<2) continue;
    int ord =  (ord1-1)*(ord2-1);
    //TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);//ponto na face
    TPZVec<FADREAL> phin(20, FADREAL(ndim, 0.0)); //3d
    //phin.Zero();
    //dphin.Zero();
    int ordin =  (ord1 > ord2) ? ord1 : ord2;
    ordin--;
    TPZManVector<int> ids(4);
    //TPZVec<int> ids(4);
    int id0,id1,i;
    for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
    id0 = ShapeFaceId[face][0];//numero das shapes da face que compoem a shape atual
    id1 = ShapeFaceId[face][1];
    TPZShapeQuad::Shape2dQuadInternal(outval,ord1-2,phin,TPZShapeQuad::GetTransformId2dQ(ids));//ordin = ordem de um lado
//    TransformDerivativeFromFaceToCube(face,ord,phin);//ord = numero de shapes
    for(i=0;i<ord;i++)	{
//      phi(shape,0) = phi(id0,0)*phi(id1,0)*phin(i,0);
        phi[shape] = phi[id0] * phi[id1] * phin[i];
/*      for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj,id0)* phi(id1 , 0 )* phin(i ,0) +
                           phi(id0, 0)*dphi(xj  ,id1)* phin(i ,0) +
                           phi(id0, 0)* phi(id1 , 0 )*dphin(xj,i);  // implicitly done
      }*/
      shape++;
    }
  }

  //volume shapes
  //REAL store1[20],store2[60];
  int ordmin1 = (order[18]-1);
  int ord =  ordmin1*ordmin1*ordmin1;//(p-1)^3 : 0<=n1,n2,n3<=p-2
  //TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);
  TPZVec<FADREAL> phin(20, FADREAL(ndim, 0.0)); //3d
  //phin.Zero();
  //dphin.Zero();
  Shape3dCubeInternal(pt,ordmin1,phin);
  for(int i=0;i<ord;i++)	{
    //phi(shape,0) = phi(0,0)*phi(6,0)*phin(i,0);
    phi[shape] = phi[0] * phi[6] * phin[i];
    /*
    for(int xj=0;xj<3;xj++) {
      dphi(xj,shape) = dphi(xj,0)* phi(6 ,0)* phin(i ,0) +
                        phi(0, 0)*dphi(xj,6)* phin(i ,0) +
                        phi(0, 0)* phi(6 ,0)*dphin(xj,i);
    }*/
    shape++;
  }
}


void TPZShapeCube::ShapeCornerCube(TPZVec<FADREAL> &pt, TPZVec<FADREAL> &phi)
{
      FADREAL x[2], y[2], z[2];

      x[0]  = (REAL(1.)-pt[0])/REAL(2.);
      x[1]  = (REAL(1.)+pt[0])/REAL(2.);
      y[0]  = (REAL(1.)-pt[1])/REAL(2.);
      y[1]  = (REAL(1.)+pt[1])/REAL(2.);
      z[0]  = (REAL(1.)-pt[2])/REAL(2.);
      z[1]  = (REAL(1.)+pt[2])/REAL(2.);

      phi[0]  = x[0]*y[0]*z[0];
      phi[1]  = x[1]*y[0]*z[0];
      phi[2]  = x[1]*y[1]*z[0];
      phi[3]  = x[0]*y[1]*z[0];
      phi[4]  = x[0]*y[0]*z[1];
      phi[5]  = x[1]*y[0]*z[1];
      phi[6]  = x[1]*y[1]*z[1];
      phi[7]  = x[0]*y[1]*z[1];
}

void TPZShapeCube::ProjectPoint3dCubeToRib(int side, TPZVec<FADREAL> &in, FADREAL &outval)
{
  outval = gRibTrans3dCube1d[side][0]*in[0]+gRibTrans3dCube1d[side][1]*in[1]+gRibTrans3dCube1d[side][2]*in[2];
/* outval =  gRibTrans3dCube1d[side][0]*in[0];
 outval += gRibTrans3dCube1d[side][1]*in[1];
 outval += gRibTrans3dCube1d[side][2]*in[2];*/
}
/*
void TPZShapeCube::TransformDerivativeFromRibToCube(int rib,int num,TPZVec<FADREAL> &phi) {
  for (int j = 0;j<num;j++) {
    //dphi(2,j) = gRibTrans3dCube1d[rib][2]*dphi(0,j);
    //dphi(1,j) = gRibTrans3dCube1d[rib][1]*dphi(0,j);
    //dphi(0,j) = gRibTrans3dCube1d[rib][0]*dphi(0,j);
    phi[j].fastAccessDx(2) = gRibTrans3dCube1d[rib][2]*phi[j].d(0);
    phi[j].fastAccessDx(1) = gRibTrans3dCube1d[rib][1]*phi[j].d(0);
    phi[j].fastAccessDx(0) = gRibTrans3dCube1d[rib][0]*phi[j].d(0);
  }
}
*/
void TPZShapeCube::ProjectPoint3dCubeToFace(int face, TPZVec<FADREAL> &in, TPZVec<FADREAL> &outval) {
  outval[0] = gFaceTrans3dCube2d[face][0][0]*in[0]+gFaceTrans3dCube2d[face][0][1]*in[1]+gFaceTrans3dCube2d[face][0][2]*in[2];
  outval[1] = gFaceTrans3dCube2d[face][1][0]*in[0]+gFaceTrans3dCube2d[face][1][1]*in[1]+gFaceTrans3dCube2d[face][1][2]*in[2];
/*  outval[0] = gFaceTrans3dCube2d[face][0][0]*in[0];
  outval[0] += gFaceTrans3dCube2d[face][0][1]*in[1];
  outval[0] += gFaceTrans3dCube2d[face][0][2]*in[2];
  outval[1] = gFaceTrans3dCube2d[face][1][0]*in[0];
  outval[1] += gFaceTrans3dCube2d[face][1][1]*in[1];
  outval[1] += gFaceTrans3dCube2d[face][1][2]*in[2];*/
}

/*
void TPZShapeCube::TransformDerivativeFromFaceToCube(int rib,int num,TPZVec<FADREAL> &phi) {

  for (int j = 0;j<num;j++) {

    //dphi(2,j) = gFaceTrans3dCube2d[rib][0][2]*dphi(0,j)+gFaceTrans3dCube2d[rib][1][2]*dphi(1,j);
    //REAL dphi1j = dphi(1,j);
    //dphi(1,j) = gFaceTrans3dCube2d[rib][0][1]*dphi(0,j)+gFaceTrans3dCube2d[rib][1][1]*dphi(1,j);
    //dphi(0,j) = gFaceTrans3dCube2d[rib][0][0]*dphi(0,j)+gFaceTrans3dCube2d[rib][1][0]*dphi1j;//dphi(1,j);

    REAL dphijd1 = phi[j].d(1);
    phi[j].fastAccessDx(2) =  gFaceTrans3dCube2d[rib][0][2]*phi[j].d(0)+gFaceTrans3dCube2d[rib][1][2]*dphijd1;
    phi[j].fastAccessDx(1) =  gFaceTrans3dCube2d[rib][0][1]*phi[j].d(0)+gFaceTrans3dCube2d[rib][1][1]*dphijd1;
    phi[j].fastAccessDx(0) =  gFaceTrans3dCube2d[rib][0][0]*phi[j].d(0)+gFaceTrans3dCube2d[rib][1][0]*dphijd1;
  }
}
*/

void TPZShapeCube::Shape3dCubeInternal(TPZVec<FADREAL> &x, int order,TPZVec<FADREAL> &phi)
{//,int cube_transformation_index
  const int ndim = 3;

  if(order < 1) return;
  int ord = order;//fSideOrder[18]-1;
  order = order*order*order;
  phi.Resize(order, FADREAL(ndim, 0.0));
  TPZVec<FADREAL> phi0(20, FADREAL(ndim, 0.0)),
                  phi1(20, FADREAL(ndim, 0.0)),
		  phi2(20, FADREAL(ndim, 0.0));
  //phi.Resize(order, 1);
  //dphi.Resize(3,order);
  //REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
  //TPZFMatrix phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
  //  dphi0(1,ord,store4,20),dphi1(1,ord,store5,20),dphi2(1,ord,store6,20);
  TPZShapeLinear::FADfOrthogonal(x[0],ord,phi0);
  TPZShapeLinear::FADfOrthogonal(x[1],ord,phi1);
  TPZShapeLinear::FADfOrthogonal(x[2],ord,phi2);
  for (int i=0;i<ord;i++) {
    for (int j=0;j<ord;j++) {
      for (int k=0;k<ord;k++) {
         int index = ord*(ord*i+j)+k;
          //phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
	  phi[index] =  phi0[i] * phi1[j] * phi2[k];
         /*dphi(0,index) = dphi0(0,i)* phi1(j,0)* phi2(k,0);
         dphi(1,index) =  phi0(i,0)*dphi1(0,j)* phi2(k,0);
         dphi(2,index) =  phi0(i,0)* phi1(j,0)*dphi2(0,k);
	 */
      }
    }
  }
}

#endif

};
