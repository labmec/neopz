// $Id: pzshapepiram.cpp,v 1.3 2003-10-06 01:32:07 phil Exp $
#include "pzshapepiram.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"


/*Projection of the point within a piramide to a rib*/
REAL TPZShapePiram::gRibTrans3dPiram1d[8][3] = {//parâmetros de arestas
  { 1., 0.,0.} , { 0., 1.,0.} ,//percorre o sentido
  {-1., 0.,0.} , { 0.,-1.,0.} ,//da aresta segundo
  { .5, .5,1.} , {-.5,.5 ,1.} ,//SideNodes[8][2]
  {-.5,-.5,1.} , { .5,-.5,1.}
};
//não tem vetor associado -> OK!

//Projection of the point within a piramide to a face
REAL TPZShapePiram::gFaceTrans3dPiram2d[5][2][3] = {//parâmetros de faces
  { { 1., 0., 0.},{ 0., 1.,0.}  },//0 1 2 3  : percorre os eixos segundo
  { { 1.,-.5,-.5},{ 0., 1.,1.}  },//0 1 4    : FaceNodes[5][4]
  { { .5, 1.,-.5},{-1., 0.,1.}  },//1 2 4
  //{ {-1., .5,-.5},{ 0.,-1.,1.}  },//ORIGINAL
  //{ {-.5,-1.,-.5},{ 1., 0.,1.}  }//ORIGINAL
  { { 1., .5,-.5},{ 0.,-1.,1.}  },//3 2 4 ; original-> 2 3 4 : {-1., .5,-.5},{ 0.,-1.,1.}
  { {-.5, 1.,-.5},{1.,0.,1.}  } //0 3 4 ; original-> 3 0 4 : {-.5,-1.,-.5},{ 1., 0.,1.}
};

int TPZShapePiram::SideNodes[8][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4} };

int TPZShapePiram::FaceNodes[5][4]  = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };

int TPZShapePiram::ShapeFaceId[5][4] = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };

REAL TPZShapePiram::gFaceSum3dPiram2d[5][2] = { {.0,0.},{-.5,0.},{-.5,0.},{-.5,0.},{-.5,0.} };//{ {.0,0.},{-.5,0.},{-.5,0.},{-.5,0.},{-.5,0.} };//original

static int sidedimension[19] = {0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,3};

static int nhighdimsides[19] = {7,7,7,7,9,3,3,3,3,3,3,3,3,1,1,1,1,1,0};

static int highsides[19][9] = {
{5,8,9,13,14,17,18},
{5,6,10,13,14,15,18},
{6,7,11,13,15,16,18},
{7,8,12,13,16,17,18},
{9,10,11,12,14,15,16,17,18},
{13,14,18},
{13,15,18},
{13,16,18},
{13,17,18},
{14,17,18},
{14,15,18},
{15,16,18},
{16,17,18},
{18},
{18},
{18},
{18},
{18},
{-999}
};

static REAL sidetosidetransforms[19][9][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
},
{
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,0}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
},
{
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
},
{
{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,0}}
},
{
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0.5,0.5,0.5},{-99,-99,-99},{-99,-99,-99},{-0.5,-0.5,0.5}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-0.5,0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,-0.5,0.5}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{-0.5,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,0.5}}
},
{
{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{0.5,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{-0.5,0.5,0.5}}
},
{
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,0}}
},
{
{{2,0,0},{1,1,1},{-99,-99,-99},{-1,-1,0}}
},
{
{{0,2,0},{-1,1,1},{-99,-99,-99},{1,-1,0}}
},
{
{{2,0,0},{1,-1,1},{-99,-99,-99},{-1,1,0}}
},
{
{{0,2,0},{1,1,1},{-99,-99,-99},{-1,-1,0}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static int FaceConnectLocId[5][9] = { {0,1,2,3,5,6,7,8,13},{0,1,4,5,10,9,14,-1,-1},
					      {1,2,4,6,11,10,15,-1,-1},{3,2,4,7,11,12,16,-1,-1},{0,3,4,8,12,9,17,-1,-1} };

static REAL MidSideNode[19][3] = {
/*00*/{-1.,-1.},   /*01*/{1.,-1.},   /*02*/{1.,1.},/*03*/{-1.,1.},/*04*/{0.,0.,1.},
/*05*/{ 0.,-1.},   /*06*/{1., 0.},   /*07*/{0.,1.},/*08*/{-1.,0.},
/*09*/{-.5,-.5,.5},/*10*/{.5,-.5,.5},/*11*/{.5,.5,.5},/*12*/{-.5,.5,.5},
/*13*/{0.,  0. ,  0. },/*14*/{  0.  ,-2./3.,1./3.},/*15*/{2./3.,0.,1./3.},
/*16*/{0.,2./3.,1./3.},/*17*/{-2./3.,  0.  ,1./3.},/*18*/{  0. ,0.,1./5.} };

void TPZShapePiram::CornerShape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {

  /*if(abs(pt[0])<1.e-10 && abs(pt[1])<1.e-10 && pt[2]==1.) {
  	//para testes com transformações geometricas
   //(0,0,1) nunca é um ponto de integração
     phi(0,0)  = 0.;
     phi(1,0)  = 0.;
     phi(2,0)  = 0.;
     phi(3,0)  = 0.;
     phi(4,0)  = 1.;
     for(int i=0;i<5;i++) {
        dphi(0,i) = 0.;
        dphi(1,i) = 0.;
        dphi(2,i) = 0.;
     }
     return;
  }*/
  REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
  REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
  REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
  REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
  REAL lmez = (1.-pt[2]);
  phi(0,0)  = T0xz*T0yz*lmez;
  phi(1,0)  = T1xz*T0yz*lmez;
  phi(2,0)  = T1xz*T1yz*lmez;
  phi(3,0)  = T0xz*T1yz*lmez;
  phi(4,0)  = pt[2];
  REAL lmexmez = 1.-pt[0]-pt[2];
  REAL lmeymez = 1.-pt[1]-pt[2];
  REAL lmaxmez = 1.+pt[0]-pt[2];
  REAL lmaymez = 1.+pt[1]-pt[2];
  dphi(0,0) = -.25*lmeymez / lmez;
  dphi(1,0) = -.25*lmexmez / lmez;
  dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;

  dphi(0,1) =  .25*lmeymez / lmez;
  dphi(1,1) = -.25*lmaxmez / lmez;
  dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;

  dphi(0,2) =  .25*lmaymez / lmez;
  dphi(1,2) =  .25*lmaxmez / lmez;
  dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;

  dphi(0,3) = -.25*lmaymez / lmez;
  dphi(1,3) =  .25*lmexmez / lmez;
  dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;

  dphi(0,4) =  0.0;
  dphi(1,4) =  0.0;
  dphi(2,4) =  1.0;
}
/*REAL TPZCompEl::T(REAL &pt0,REAL &pt1) {
  return ( ( ( 1 - pt1 ) - pt0 ) / 2 / (1 - pt1) );
}

void TPZCompEl::ShapeCornerPira(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {


  double c = pt[0] , e = pt[1], z = 1 - pt[2];
  phi(0,0) = T(c,pt[2]) * T( e,pt[2]) * z;
  c *= -1;
  phi(1,0) = T(c,pt[2]) * T( e,pt[2]) * z;
  e *= -1;
  phi(2,0) = T(c,pt[2]) * T(e,pt[2]) * z;
  c *= -1;
  phi(3,0) = T(c,pt[2]) * T(e,pt[2]) * z;
  phi(4,0) = pt[2];
  e *= -1;

  dphi(0,0) = -.25 * (z - e)/z;
  dphi(1,0) = -.25 * (z - c)/z;
  dphi(2,0) = -.25 * (z - c)/z - .25 * (z - e)/z + .25 * (z - c) * (z - e) / z / z ;

  dphi(0,1) =  .25 * (z - e)/z;
  dphi(1,1) = -.25 * (z + c)/z;
  dphi(2,1) = -.25 * (z + c)/z - .25 * (z - e)/z + .25 * (z + c) * (z - e) / z / z ;

  dphi(0,2) =  .25 * (z + e)/z;
  dphi(1,2) =  .25 * (z + c)/z;
  dphi(2,2) = -.25 * (z + c)/z - .25 * (z + e)/z + .25 * (z + c) * (z + e) / z / z ;

  dphi(0,3) = -.25 * (z + e)/z;
  dphi(1,3) =  .25 * (z - c)/z;
  dphi(2,3) = -.25 * (z - c)/z - .25 * (z + e)/z + .25 * (z - c) * (z + e) / z / z ;

  dphi(0,4) = 0.0;
  dphi(1,4) = 0.0;
  dphi(2,4) = 1.0;
}*/

void TPZShapePiram::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {

  CornerShape(pt,phi,dphi);
  //  if(order[13]<2) return;//order tem as ordens dos lados do elemento
  int shape = 5;
  //rib shapes
  for (int rib = 0; rib < 8; rib++) {//todas as arestas
    if (order[rib] <2 ) continue;
    REAL outval;
    ProjectPoint3dPiramToRib(rib,pt,outval);
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
    TransformDerivativeFromRibToPiram(rib,ordin,dphin);
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
  //if(order[13]<2) return;//ordem do elemento
  //face shapes
  for (int face = 0; face < 5; face++) {
    if (order[face+8] < 2) continue;
    if(face>0 && order[face+8]<=2) continue;//só a face 13 tem shape associada com ordem p=2
    TPZVec<REAL> outval(2);
    ProjectPoint3dPiramToFace(face,pt,outval);
    REAL store1[20],store2[60];
    int ord1;//,ord2;
	 ord1 = order[8+face];
	//ord2 = ord1;
    //FaceOrder(face,ord1,ord2);//ordem da face
    if(face && ord1<3) continue;//uma face com ordem < 3 não tem shape associada
    int ordin;
    if(!face) ordin = (ord1-1)*(ord1-1);//face quadrada
    else ordin = (ord1-2)*(ord1-1)/2;//face triangular
    TPZFMatrix phin(ordin,1,store1,20),dphin(3,ordin,store2,60);//ponto na face
    phin.Zero();
    dphin.Zero();
    TPZManVector<int> ids(4);
    int id0,id1,id2,i;
    if(!face) for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
	 else for(i=0;i<3;i++) ids[i] = id[FaceNodes[face][i]];
    id0 = ShapeFaceId[face][0];//indice das shapes da face que compoem a shape atual
    id1 = ShapeFaceId[face][1];//equivale a FaceIdsCube(face,ids,id,id0,id1);
    id2 = ShapeFaceId[face][2];//if(face == 0) id3 = ShapeFaceId[face][3];
    int transid;
    if(!face) {
      transid = TPZShapeQuad::GetTransformId2dQ(ids);
      TPZShapeQuad::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
    } else {
       ids.Resize(3);
    	 transid = TPZShapeTriang::GetTransformId2dT(ids);
       outval[0] = (outval[0]+1.)/2.;//devido a correção na função
       outval[1] = (outval[1]+1.)/2.;//Shape2dTriangleInternal(..) : correto aqui
       TPZShapeTriang::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
       int c = dphin.Cols();//isto da (p-2)(p-1)/2 ; ord1 = p ; correto aqui
       for(int i=0;i<c;i++) {
         dphin(0,i) /= 2.;//correcao da derivada OK! aqui
         dphin(1,i) /= 2.;
         dphin(2,i) /= 2.;
       }
    }
    TransformDerivativeFromFaceToPiram(face,ordin,dphin);//ordin = numero de shapes
    for(i=0;i<ordin;i++)	{
      phi(shape,0) = phi(id0,0)*phi(id2,0)*phin(i,0);//face quadriláteral
      REAL fi1 = phi(id1 , 0 );
      if(face) phi(shape,0) *= fi1;//face triangular
      for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj,id0)* phi(id2 , 0 )* phin(i ,0) +
                           phi(id0, 0)*dphi(xj  ,id2)* phin(i ,0) +
                           phi(id0, 0)* phi(id2 , 0 )*dphin(xj,i);

         if(face) {
         	dphi(xj,shape) *= fi1;
            dphi(xj,shape) += phi(id0, 0)* phi(id2 , 0 )*dphi(xj  ,id1)* phin(i ,0);
         }
      }
      shape++;
    }
  }
  if(order[13]<3) return;//não há ordens para cantos
  //volume shapes
  REAL store1[20],store2[60];
  int ord=0,i;
  for(i=0;i<order[13]-2;i++) {
    ord += (i+1)*(i+2) / 2;
  }
  TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);
  phin.Zero();
  dphin.Zero();
  ShapeInternal(pt,order[13],phin,dphin);
  for(i=0;i<ord;i++)	{
    phi(shape,0) = phi(0,0)*phi(2,0)*phi(4,0)*phin(i,0);
    for(int xj=0;xj<3;xj++) {
      dphi(xj,shape) = dphi(xj,0)* phi(2 ,0)* phi(4 ,0)* phin(i ,0) +
                        phi(0, 0)*dphi(xj,2)* phi(4 ,0)* phin(i ,0) +
                        phi(0, 0)* phi(2 ,0)*dphi(xj,4)* phin(i ,0) +
                        phi(0, 0)* phi(2 ,0)* phi(4 ,0)* dphin(xj,i);
    }
    shape++;
  }
}

void TPZShapePiram::SideShape(int side, TPZVec<REAL> &point, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {

  if(side<0 || side>18) PZError << "TPZCompElPi3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==18) Shape(point,id,order,phi,dphi);
  else if(side<5) TPZShapePoint::Shape(point,id,order,phi,dphi);
  else if(side<13) {//5 a 12
    TPZShapeLinear::Shape(point,id,order,phi,dphi);
  } else if(side == 13) {
    TPZShapeQuad::Shape(point,id,order,phi,dphi);
  } else if(side<18) {//faces 13  a 17
      TPZShapeTriang::Shape(point,id,order,phi,dphi);
  }

}

void TPZShapePiram::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                    TPZFMatrix &dphi) {
  //valor da função e derivada do produto das funções ortogonais
  if(order < 3) return;
  int ord = order-2;
  REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
  TPZFMatrix phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
            dphi0(1,ord,store4,20),dphi1(1,ord,store5,20),dphi2(1,ord,store6,20);
  TPZShapeLinear::fOrthogonal(x[0],ord,phi0,dphi0);//f e df            -1<=x0<=1
  TPZShapeLinear::fOrthogonal(x[1],ord,phi1,dphi1);//g e dg            -1<=x1<=1
  TPZShapeLinear::fOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);//h e dh       0<=x3<=1 -> -1<=2*x2-1<=1
  int index = 0;//x é ponto de integraçào dentro da pirâmide
  for (int i=0;i<ord;i++) {
    for (int j=0;j<ord;j++) {
      for (int k=0;k<ord;k++) {
      	if( i+j+k < ord ) {
            //int index = ord*(ord*i+j)+k; //i,j,k é o grau das funções ortogonais
             phi(index,0) =    phi0(i,0)* phi1(j,0)* phi2(k,0);
            dphi(0,index) =   dphi0(0,i)* phi1(j,0)* phi2(k,0);
            dphi(1,index) =    phi0(i,0)*dphi1(0,j)* phi2(k,0);
            dphi(2,index) = 2.*phi0(i,0)* phi1(j,0)*dphi2(0,k);
            index++;
         }
      }
    }
  }
}

void TPZShapePiram::TransformDerivativeFromRibToPiram(int rib,int num,TPZFMatrix &dphi) {
  for (int j = 0;j<num;j++) {
    dphi(2,j) = gRibTrans3dPiram1d[rib][2]*dphi(0,j);
    dphi(1,j) = gRibTrans3dPiram1d[rib][1]*dphi(0,j);
    dphi(0,j) = gRibTrans3dPiram1d[rib][0]*dphi(0,j);
  }
}

void TPZShapePiram::TransformDerivativeFromFaceToPiram(int face,int num,TPZFMatrix &dphi) {
  for (int j = 0;j<num;j++) {
    dphi(2,j) = gFaceTrans3dPiram2d[face][0][2]*dphi(0,j)+gFaceTrans3dPiram2d[face][1][2]*dphi(1,j);
    REAL dphi1j = dphi(1,j);
    dphi(1,j) = gFaceTrans3dPiram2d[face][0][1]*dphi(0,j)+gFaceTrans3dPiram2d[face][1][1]*dphi(1,j);
    dphi(0,j) = gFaceTrans3dPiram2d[face][0][0]*dphi(0,j)+gFaceTrans3dPiram2d[face][1][0]*dphi1j;//dphi(1,j);
  }
}
//transforma a derivada no ponto dentro da face
void TPZShapePiram::TransformDerivativeFace3dPiram(int transid, int face, int num, TPZFMatrix &in) {
  if (!face) TPZShapeQuad::TransformDerivative2dQ(transid,num,in);//face 13
  else       TPZShapeTriang::TransformDerivative2dT(transid,num,in);//outras
}

//projeta o ponto do interior para o lado
void TPZShapePiram::ProjectPoint3dPiramToRib(int rib, TPZVec<REAL> &in, REAL &outval) {
  outval = gRibTrans3dPiram1d[rib][0]*in[0]+gRibTrans3dPiram1d[rib][1]*in[1]+gRibTrans3dPiram1d[rib][2]*in[2];
}

//projeta o ponto do interior para a face
void TPZShapePiram::ProjectPoint3dPiramToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval) {
  outval[0] = gFaceTrans3dPiram2d[face][0][0]*in[0]+gFaceTrans3dPiram2d[face][0][1]*in[1]+gFaceTrans3dPiram2d[face][0][2]*in[2]+gFaceSum3dPiram2d[face][0];
  outval[1] = gFaceTrans3dPiram2d[face][1][0]*in[0]+gFaceTrans3dPiram2d[face][1][1]*in[1]+gFaceTrans3dPiram2d[face][1][2]*in[2]+gFaceSum3dPiram2d[face][1];
}

//transforma o ponto dentro da face
void TPZShapePiram::TransformPoint3dPiramFace(int transid, int face, TPZVec<REAL> &in, TPZVec<REAL> &out) {
  if (!face) TPZShapeQuad::TransformPoint2dQ(transid,in,out);//face zero ou 13
  else       TPZShapeTriang::TransformPoint2dT(transid,in,out);//outras 14 a 17
}

int TPZShapePiram::NConnectShapeF(int side, int order) {
   if(side<5) return 1;//0 a 4
   //   int s = side-5;//s = 0 a 14 ou side = 6 a 20
   if(side<13) return (order-1);//6 a 14
   if(side==13) {
      return ((order-1)*(order-1));
   }
   if(side<18) {//16,17,18
      return ((order-2)*(order-1)/2);
   }
   if(side==18) {
   	int totsum = 0,sum;
      for(int i=1;i<order-1;i++) {
         sum = i*(i+1) / 2;
         totsum += sum;
      }
      return totsum;
   }
   PZError << "TPZShapePiram::NConnectShapeF, bad parameter side " << side << endl;
   return 0;
}

int TPZShapePiram::NShapeF(TPZVec<int> &order) {
  int in,res=NNodes;
  for(in=NNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NNodes]);
  return res;
}

int TPZShapePiram::NConnects() {
	return 19;
}

TPZTransform TPZShapePiram::TransformElementToSide(int side){

	if(side<0 || side>18){
  	PZError << "TPZShapePiram::TransformElementToSide called with side error\n";
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
      return t;
    case 5:
    	t.Mult()(0,0) = 1.0;
      return t;
    case 6:
      t.Mult()(0,1) = 1.0;
      return t;
    case 7:
    	t.Mult()(0,0) = -1.0;
      return t;
    case 8:
    	t.Mult()(0,1) = -1.0;
      return t;
    case 9:
    case 12:
    	t.Mult()(0,0) = 2.0;
      t.Sum()(0,0)  = 1.0;
      return t;
    case 10:
    case 11:
    	t.Mult()(0,0) = -2.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case 13:
    	t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
    case 14:
    	t.Mult()(0,0) =  0.5;
      t.Mult()(0,1) = -0.5;
      t.Mult()(1,2) =  1.0;
      return t;
    case 15:
    case 16:
    	t.Mult()(0,0) =  0.5;
      t.Mult()(0,1) =  0.5;
      t.Mult()(1,2) =  1.0;
      return t;
    case 17:
      t.Mult()(0,0) = -0.5;
    	t.Mult()(0,1) =  0.5;
      t.Mult()(1,2) =  1.0;
      return t;
    case 18:
    	t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapePiram::TransformSideToElement(int side){

	if(side<0 || side>18){
  	PZError << "TPZShapePiram::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(3,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) = -1.0;
      t.Sum()(2,0) =  0.0;
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) = -1.0;
      t.Sum()(2,0) =  0.0;
      return t;
    case 2:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) =  0.0;
      return t;
    case 3:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) =  0.0;
      return t;
    case 4:
    	t.Sum()(0,0) =  0.0;
      t.Sum()(1,0) =  0.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 5:
      t.Mult()(0,0) =  1.0;
      t.Sum() (1,0) = -1.0;
      return t;
    case 6:
      t.Mult()(1,0) =  1.0;
      t.Sum() (0,0) =  1.0;
      return t;
    case 7:
      t.Mult()(0,0) = -1.0;
      t.Sum() (1,0) =  1.0;
      return t;
    case 8:
      t.Mult()(1,0) = -1.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 9:
    	t.Mult()(0,0) =  0.5;
      t.Mult()(1,0) =  0.5;
      t.Mult()(2,0) =  0.5;
    	t.Sum()(0,0)  = -0.5;
      t.Sum()(1,0)  = -0.5;
      t.Sum()(2,0)  =  0.5;
      return t;
    case 10:
    	t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Mult()(2,0) =  0.5;
    	t.Sum()(0,0)  =  0.5;
      t.Sum()(1,0)  = -0.5;
      t.Sum()(2,0)  =  0.5;
      return t;
    case 11:
    	t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) = -0.5;
      t.Mult()(2,0) =  0.5;
    	t.Sum()(0,0)  =  0.5;
      t.Sum()(1,0)  =  0.5;
      t.Sum()(2,0)  =  0.5;
      return t;
    case 12:
    	t.Mult()(0,0) =  0.5;
      t.Mult()(1,0) = -0.5;
      t.Mult()(2,0) =  0.5;
    	t.Sum()(0,0)  = -0.5;
      t.Sum()(1,0)  =  0.5;
      t.Sum()(2,0)  =  0.5;
      return t;
    case 13:
    	t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
    case 14:
			t.Mult()(0,0) =  2.0;
    	t.Mult()(0,1) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,1) =  1.0;
    	t.Sum()(0,0)  = -1.0;
      t.Sum()(1,0)  = -1.0;
      return t;
    case 15:
			t.Mult()(1,0) =  2.0;
    	t.Mult()(0,1) = -1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,1) =  1.0;
    	t.Sum()(0,0)  =  1.0;
      t.Sum()(1,0)  = -1.0;
      return t;
    case 16:
			t.Mult()(0,0) =  2.0;
    	t.Mult()(0,1) =  1.0;
      t.Mult()(1,1) = -1.0;
      t.Mult()(2,1) =  1.0;
    	t.Sum()(0,0)  = -1.0;
      t.Sum()(1,0)  =  1.0;
      return t;
    case 17:
			t.Mult()(1,0) =  2.0;
    	t.Mult()(0,1) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,1) =  1.0;
    	t.Sum()(0,0)  = -1.0;
      t.Sum()(1,0)  = -1.0;
      return t;
    case 18:
			t.Mult()(0,0) =  1.0;
    	t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}


int TPZShapePiram::SideNodeLocId(int side, int node)
{
	if(side <5 && node == 0) return side;
	if(side >= 5 && side <13 && node < 2) return SideNodes[side-5][node];
	if(side == 13 && node <4) return FaceNodes[side-13][node];
	if(side >13 && side < 18)
	  if (node <3) return FaceNodes[side-13][node];
	  else if (node==3) return -1;//Previsto receber pelas faces triangulares - Cesar 2003-01-02
	
	if(side == 18 && node < 5) return node;
	PZError << "TPZShapePiram::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

int nsidenodes[19] = {1,1,1,1,1,
2,2,2,2,2,2,2,2,
4,3,3,3,3,
5};

int TPZShapePiram::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapePiram::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapePiram::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapePiram::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapePiram::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZShapePiram::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapePiram::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapePiram::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

int TPZShapePiram::NSideConnects(int side) {
  if(side<0)   return -1;
  if(side<5)   return 1;//cantos : 0 a 4
  if(side<13)  return 3;//lados : 5 a 12
  if(side==13) return 9;//face : 13 , quadrilateral
  if(side<18)	 return 7;//faces : 14 a 17 , triangulares
  if(side==18) return 19;//centro : 18
  return -1;
}

int TPZShapePiram::SideConnectLocId(int side, int node) {
  if(side<0 || side>19 || node < 0) return -1;
  if(side<5) {
    if(node==0) return side;
  } 
  else if(side<9) {//5 a 8
    int s = side-5;//0,1,2
    if(!node) return s;//0,1,2,3
    if(node==1) return (s+1)%4;//1,2,0
    if(node==2) return side;//5,6,7,8
  }
  else if(side<13) {//9 a 12
    int s = side-9;//0,1,2,3
    if(!node) return s;//0,1,2,3
    if(node==1) return 4;//
    if(node==2) return side;//9,10,11,12
  } 
  else if(side<18) {//13 a 17
    int s = side-13;
    if(node<9) return FaceConnectLocId[s][node];
  } 
  else if(side==18 && node<19){
    return node;
  }
  PZError << "TPZShapePiram::SideConnectLocId called for node = "
	  << node << " and side = " << side << "\n";
  return -1;
}

void TPZShapePiram::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}
