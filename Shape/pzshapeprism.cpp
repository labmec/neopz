// $Id: pzshapeprism.cpp,v 1.4 2003-11-25 17:59:45 cesar Exp $
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"


/*Projection of the point within a piramide to a rib*/
REAL TPZShapePrism::gRibTrans3dPrisma1d[9][3] = {//parâmetros de arestas
  { 2., 1.,0.} , {-1., 1.,0.} ,//percorre o sentido
  {-1.,-2.,0.} , { 0., 0.,1.} ,//da aresta segundo : { 1., 2.,0.} , { 0., 0.,1.}
  { 0., 0.,1.} , { 0., 0.,1.} ,//SideNodes[9][2]
  { 2., 1.,0.} , {-1., 1.,0.} ,
  {-1.,-2.,0.}                                     //{ 1., 2.,0.}
};

REAL TPZShapePrism::gRibSum3dPrisma1d[9] = {-1.,0.,1.,0.,0.,0.,-1.,0.,1.};;//{-1.,0.,-1.,0.,0.,0.,-1.,0.,-1.};

//Projection of the point within a piramide to a face
REAL TPZShapePrism::gFaceTrans3dPrisma2d[5][2][3] = {//parâmetros de faces
  { { 2., 0., 0.},{ 0., 2., 0.}  },//0 1 2   : percorre os eixos segundo
  { { 2., 0., 0.},{ 0., 0., 1.}  },//0 1 4 3    : FaceNodes[5][4]
  { {-1., 1., 0.},{ 0., 0., 1.}  },//1 2 5 4
  { { 0., 2., 0.},{ 0., 0., 1.}  },//0 2 5 3
  { { 2., 0., 0.},{ 0., 2., 0.}  } //3 4 5
};

REAL TPZShapePrism::gFaceSum3dPrisma2d[5][2] = { {-1.,-1.},{-1.,0.},{0.,0.},{-1.,0.},{-1.,-1.} };

int TPZShapePrism::FaceNodes[5][4]  = { {0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1} };
                                         //F15        F16       F17       F18        F19

int TPZShapePrism::SideNodes[9][2]  = { {0,1},{1,2},{2,0},{0,3},{1,4},{2,5},{3,4},{4,5},{5,3} };
                              //arestas   6     7      8    9     10    11    12    13    14

int TPZShapePrism::ShapeFaceId[5][4] = { {0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1} };
                                          //F15        F16       F17       F18       F19

static int sidedimension[21] = {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3};

static int nhighdimsides[21] = {7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,1,1,1,1,1,0};

static int highsides[21][7] = {
{6,8,9,15,16,18,20},
{6,7,10,15,16,17,20},
{7,8,11,15,17,18,20},
{9,12,14,16,18,19,20},
{10,12,13,16,17,19,20},
{11,13,14,17,18,19,20},
{15,16,20},
{15,17,20},
{15,18,20},
{16,18,20},
{16,17,20},
{17,18,20},
{16,19,20},
{17,19,20},
{18,19,20},
{20},
{20},
{20},
{20},
{20},
{-999}
};

static REAL sidetosidetransforms[21][7][4][3] = {
{
//0
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-1}}
},
{
//1
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
//{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},//está errado deve ser {-1,-1,-99}
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},//CEDRIC
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-1}}
},
{
//2
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-1}}
},
{
//3
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
},
{
//4
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,1}}
},
{
//5
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,1}}
},
{
//6
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,-1}}
},
{
//7
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-1}}
},
{
//8
{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,-1}}
},
{
//9
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{0,0,0}}
},
{
//10
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
},
{
//11
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
{{0,0,1},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
},
{
//12
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,1}}
},
{
//13
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,1}}
},
{
//14
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,1}}
},
{
//15
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,-1}}
},
{
//16
{{0.5,0,0},{0,0,1},{-99,-99,-99},{0.5,0,0}}
},
{
//17
{{-0.5,0.5,0},{0,0,1},{-99,-99,-99},{0.5,0.5,0}}
},
{
//18
{{0,0.5,0},{0,0,1},{-99,-99,-99},{0,0.5,0}}
},
{
//19
{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,1}}
},
{
//20
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static int FaceConnectLocId[5][9] = { {0,1,2,6,7,8,15,-1,-1},{0,1,4,3,6,10,12,9,16},
					      {1,2,5,4,7,11,13,10,17},{0,2,5,3,8,11,14,9,18},{3,4,5,12,13,14,19,-1,-1} };

static REAL MidSideNode[21][3] = {
/*00*/{0.,.0,-1.},/*01*/{1.,0.,-1.},/*02*/{.0,1.,-1.},/*03*/{.0,.0, 1.},
/*04*/{1.,.0, 1.},/*05*/{0.,1., 1.},/*06*/{.5,.0,-1.},/*07*/{.5,.5,-1.},
/*08*/{.0,.5,-1.},/*09*/{0.,.0, 0.},/*10*/{1.,.0, 0.},/*11*/{.0,1., 0.},
/*12*/{.5,.0, 1.},/*13*/{.5,.5, 1.},/*14*/{.0,.5, 1.},/*15*/{1./3.,1./3.,-1.},
/*16*/{.5,.0, 0.},/*17*/{.5,.5, 0.},/*18*/{0.,.5, 0.},/*19*/{1./3.,1./3., 1.},
/*20*/{1./3.,1./3.,0.} };

void TPZShapePrism::CornerShape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  phi(0,0)  = .5*(1.-pt[0]-pt[1])*(1.-pt[2]);
  phi(1,0)  = .5*pt[0]*(1.-pt[2]);
  phi(2,0)  = .5*pt[1]*(1.-pt[2]);
  phi(3,0)  = .5*(1.-pt[0]-pt[1])*(1.+pt[2]);
  phi(4,0)  = .5*pt[0]*(1.+pt[2]);
  phi(5,0)  = .5*pt[1]*(1.+pt[2]);

  dphi(0,0) = -.5*(1.-pt[2]);
  dphi(1,0) = -.5*(1.-pt[2]);
  dphi(2,0) = -.5*(1.-pt[0]-pt[1]);

  dphi(0,1) =  .5*(1.-pt[2]);
  dphi(1,1) =  .0;
  dphi(2,1) = -.5*pt[0];

  dphi(0,2) =  .0;
  dphi(1,2) =  .5*(1.-pt[2]);
  dphi(2,2) = -.5*pt[1];

  dphi(0,3) = -.5*(1.+pt[2]);
  dphi(1,3) = -.5*(1.+pt[2]);
  dphi(2,3) =  .5*(1.-pt[0]-pt[1]);

  dphi(0,4) =  .5*(1.+pt[2]);
  dphi(1,4) =  .0;
  dphi(2,4) =  .5*pt[0];

  dphi(0,5) =  .0;
  dphi(1,5) =  .5*(1.+pt[2]);
  dphi(2,5) =  .5*pt[1];
}

void TPZShapePrism::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {

  CornerShape(pt,phi,dphi);
  //  if(order[14]<2) return;//order tem as ordens dos lados do elemento
  int shape = 6;
  //rib shapes
  for (int rib = 0; rib < 9; rib++) {//todas as arestas
    REAL outval;
    ProjectPoint3dPrismaToRib(rib,pt,outval);
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
    TransformDerivativeFromRibToPrisma(rib,ordin,dphin);
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
  //  if(order[14]<2) return;//ordem do elemento
  //face shapes
  for (int face = 0; face < 5; face++) {

    if((face==0 || face==4) && order[face+9]==2) continue;//estas face nao tem shapes associadas com ordem p=2
    TPZVec<REAL> outval(2);
    ProjectPoint3dPrismaToFace(face,pt,outval);
    REAL store1[20],store2[60];
    int ord1;//,ord2;
	ord1 = order[face+9];
	//ord2 = ord1;
    //elpr->FaceOrder(face,ord1,ord2);//ordem da face
    if((face==0 || face==4) && ord1<3) continue;//uma face triangular com ordem < 3 não tem shape associada
    int ordin;
    if(face && face<4) ordin = (ord1-1)*(ord1-1);//faces quadrilaterais
    else ordin = (ord1-2)*(ord1-1)/2;//face triangular
    TPZFMatrix phin(ordin,1,store1,20),dphin(3,ordin,store2,60);//ponto na face
    phin.Zero();
    dphin.Zero();
    TPZManVector<int> ids(4);
    int id0,id1,id2,i;
    if(face ==0 || face == 4) for(i=0;i<3;i++) ids[i] = id[FaceNodes[face][i]];
	else for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
    id0 = ShapeFaceId[face][0];//indice das shapes da face x
    id1 = ShapeFaceId[face][1];//que compoem a shape atual
    id2 = ShapeFaceId[face][2];
    int transid;
    if(face && face<4) {
      transid = TPZShapeQuad::GetTransformId2dQ(ids);
    	TPZShapeQuad::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
    } else {
       ids.Resize(3);
    	 transid = TPZShapeTriang::GetTransformId2dT(ids);
       outval[0] = (outval[0]+1.)/2.;//devido a correção na função
       outval[1] = (outval[1]+1.)/2.;//Shape2dTriangleInternal(..) : correto aqui
    	 TPZShapeTriang::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
       int c = dphin.Cols();//isto da (p-2)(p-1)/2 ; ord1 = p ; correto aqui
       for(i=0;i<c;i++) {
         dphin(0,i) /= 2.;//correcao da derivada OK! aqui
         dphin(1,i) /= 2.;
         //dphin(2,i) /= 2.;
       }
    }
    TransformDerivativeFromFaceToPrisma(face,ordin,dphin);//ordin = numero de shapes
    for(i=0;i<ordin;i++)	{
      phi(shape,0) = phi(id0,0)*phi(id2,0)*phin(i,0);//face quadriláteral
      REAL fi1 = phi(id1 , 0 );
      if(face==0 || face==4) phi(shape,0) *= fi1;//face triangular
      for(int xj=0;xj<3;xj++) {
         dphi(xj,shape) = dphi(xj,id0)* phi(id2 , 0 )* phin(i ,0) +
                           phi(id0, 0)*dphi(xj  ,id2)* phin(i ,0) +
                           phi(id0, 0)* phi(id2 , 0 )*dphin(xj,i);

         if(face==0 || face==4) {
         	dphi(xj,shape) *= fi1;
            dphi(xj,shape) += phi(id0, 0)* phi(id2 , 0 )*dphi(xj  ,id1)* phin(i ,0);
         }
      }
      shape++;
    }
  }
  if(order[14]<3) return;
  //volume shapes
  REAL store1[20],store2[60];
  int ord=0;
  ord = NConnectShapeF(20,order[20-NNodes]);
  TPZFMatrix phin(ord,1,store1,20),dphin(3,ord,store2,60);
  phin.Zero();
  dphin.Zero();
  ShapeInternal(pt,order[14],phin,dphin);
  for(int i=0;i<ord;i++)	{
    phi(shape,0) = phi(0,0)*phi(4,0)*pt[1]*phin(i,0);
    for(int xj=0;xj<3;xj++) {
      dphi(xj,shape) = dphi(xj,0)* phi(4 ,0)*   pt[1]  * phin(i ,0) +
                        phi(0, 0)*dphi(xj,4)*   pt[1]  * phin(i ,0) +
                        phi(0, 0)* phi(4 ,0)*   pt[1]  * dphin(xj,i);

      if(xj == 1) dphi(xj,shape) += phi(0, 0)* phi(4 ,0)* 1. * phin(i ,0);
    }
    shape++;
  }

//}//while

}

void TPZShapePrism::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi) {
  if(side<0 || side>20) PZError << "TPZCompElPr3d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==20) {
    Shape(pt,id,order,phi,dphi);
  } else if(side<6) {
    TPZShapePoint::Shape(pt,id,order,phi,dphi);
  } else if(side < 15) {
    TPZShapeLinear::Shape(pt,id,order,phi,dphi);
  } else if(side == 15 || side == 19) {
    TPZShapeTriang::Shape(pt,id,order,phi,dphi);
  } else {
    TPZShapeQuad::Shape(pt,id,order,phi,dphi);
  }
}


void TPZShapePrism::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                    TPZFMatrix &dphi) {
  //valor da função e derivada do produto das funções ortogonais
  if(order < 3) return;
  int ord1 = order-2;
  int ord2 = order-1;
  REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
  TPZFMatrix phi0(ord1,1,store1,20),phi1(ord1,1,store2,20),phi2(ord2,1,store3,20),
            dphi0(1,ord1,store4,20),dphi1(1,ord1,store5,20),dphi2(1,ord2,store6,20);
  TPZShapeLinear::fOrthogonal(2.*x[0]-1.,ord1,phi0,dphi0);//f e df       0<=x0<=1 -> -1<=2*x0-1<=1
  TPZShapeLinear::fOrthogonal(2.*x[1]-1.,ord1,phi1,dphi1);//g e dg             0<=x1<=1 -> -1<=2*x1-1<=1
  TPZShapeLinear::fOrthogonal(x[2],ord2,phi2,dphi2);//h e dh      -1<=x3<=1
  int index = 0;//x é ponto de integraçào dentro da pirâmide
  for (int i=0;i<ord1;i++) {
    for (int j=0;j<ord1;j++) {
      for (int k=0;k<ord2;k++) {
      	if( i+j < ord1 && k < ord2) {
             phi(index,0) =     phi0(i,0)* phi1(j,0)* phi2(k,0);
            dphi(0,index) = 2.*dphi0(0,i)* phi1(j,0)* phi2(k,0);
            dphi(1,index) =  2.*phi0(i,0)*dphi1(0,j)* phi2(k,0);
            dphi(2,index) =     phi0(i,0)* phi1(j,0)*dphi2(0,k);
            index++;
         }
      }
    }
  }
}

void TPZShapePrism::TransformDerivativeFromRibToPrisma(int rib,int num,TPZFMatrix &dphi) {
  for (int j = 0;j<num;j++) {
    dphi(2,j) = gRibTrans3dPrisma1d[rib][2]*dphi(0,j);
    dphi(1,j) = gRibTrans3dPrisma1d[rib][1]*dphi(0,j);
    dphi(0,j) = gRibTrans3dPrisma1d[rib][0]*dphi(0,j);
  }
}

void TPZShapePrism::TransformDerivativeFromFaceToPrisma(int face,int num,TPZFMatrix &dphi) {
  for (int j = 0;j<num;j++) {
    dphi(2,j) = gFaceTrans3dPrisma2d[face][0][2]*dphi(0,j)+gFaceTrans3dPrisma2d[face][1][2]*dphi(1,j);
    REAL dphi1j = dphi(1,j);
    dphi(1,j) = gFaceTrans3dPrisma2d[face][0][1]*dphi(0,j)+gFaceTrans3dPrisma2d[face][1][1]*dphi(1,j);
    dphi(0,j) = gFaceTrans3dPrisma2d[face][0][0]*dphi(0,j)+gFaceTrans3dPrisma2d[face][1][0]*dphi1j;//dphi(1,j);
  }
}
//transforma a derivada no ponto dentro da face
void TPZShapePrism::TransformDerivativeFace3dPrisma(int transid, int face, int num, TPZFMatrix &in) {
  if (face==0 || face==4) TPZShapeTriang::TransformDerivative2dT(transid,num,in);
  else                    TPZShapeQuad::TransformDerivative2dQ(transid,num,in);

}

//projeta o ponto do interior para o lado
void TPZShapePrism::ProjectPoint3dPrismaToRib(int rib, TPZVec<REAL> &in, REAL &outval) {
  outval = gRibTrans3dPrisma1d[rib][0]*in[0]+gRibTrans3dPrisma1d[rib][1]*in[1]+gRibTrans3dPrisma1d[rib][2]*in[2]+gRibSum3dPrisma1d[rib];
}

//projeta o ponto do interior para a face
void TPZShapePrism::ProjectPoint3dPrismaToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval) {
  outval[0] = gFaceTrans3dPrisma2d[face][0][0]*in[0]+gFaceTrans3dPrisma2d[face][0][1]*in[1]+gFaceTrans3dPrisma2d[face][0][2]*in[2]+gFaceSum3dPrisma2d[face][0];
  outval[1] = gFaceTrans3dPrisma2d[face][1][0]*in[0]+gFaceTrans3dPrisma2d[face][1][1]*in[1]+gFaceTrans3dPrisma2d[face][1][2]*in[2]+gFaceSum3dPrisma2d[face][1];
}

//transforma o ponto dentro da face
void TPZShapePrism::TransformPoint3dPrismaFace(int transid, int face, TPZVec<REAL> &in, TPZVec<REAL> &out) {
  if (face==0 || face==4) TPZShapeTriang::TransformPoint2dT(transid,in,out);
  else                    TPZShapeQuad::TransformPoint2dQ(transid,in,out);
}

int TPZShapePrism::NConnectShapeF(int side, int order) {
   if(side<6) return 1;//0 a 4
   //   int s = side-6;//s = 0 a 14 ou side = 6 a 20
   if(side<15) return (order-1);//6 a 14
   if(side==15 || side==19) {
      return ((order-2)*(order-1)/2);
   }
   if(side>15 && side<19) {//16,17,18
      return ((order-1)*(order-1));
   }
   if(side==20) {
      return ((order-2)*(order-1)*(order-1)/2);
   }
   PZError << "TPZShapePrism::NConnectShapeF, bad parameter side " << side << endl;
   return 0;
}

int TPZShapePrism::NConnects() {
	return 21;
}

int TPZShapePrism::NShapeF(TPZVec<int> &order) {
  int in,res=NNodes;
  for(in=NNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NNodes]);
  return res;
}

TPZTransform TPZShapePrism::TransformElementToSide(int side){

	if(side<0 || side>20){
  	PZError << "TPZShapePrism::TransformElementToSide called with side error\n";
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
    case 5:
      return t;
    case  6:
    case 12:
    	t.Mult()(0,0) =  2.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case  7:
    case 13:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) =  1.0;
      return t;
    case  8:
    case 14:
    	t.Mult()(0,1) = -2.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case  9:
    case 10:
    case 11:
    	t.Mult()(0,2) =  1.0;
      return t;
    case 15:
    case 19:
    	t.Mult()(0,0) = 1.0;
      t.Mult()(1,1) = 1.0;
      return t;
    case 16:
    	t.Mult()(0,0) =  2.0;
      t.Mult()(1,2) =  1.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 17:
    	t.Mult()(0,0) = -1.0;
        t.Mult()(0,1) =  1.0;
	    t.Mult()(1,2) =  1.0;
      return t;
    case 18:
    	t.Mult()(0,1) =  2.0;
      t.Mult()(1,2) =  1.0;
      t.Sum()(0,0)  = -1.0;
      return t;
    case 20:
    	t.Mult()(0,0) = 1.0;
      t.Mult()(1,1) = 1.0;
      t.Mult()(2,2) = 1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapePrism::TransformSideToElement(int side){

	if(side<0 || side>20){
  	PZError << "TPZShapePrism::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(3,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
    	t.Sum()(0,0) =  0.0;
      t.Sum()(1,0) =  0.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  0.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 2:
    	t.Sum()(0,0) =  0.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) = -1.0;
      return t;
    case 3:
    	t.Sum()(0,0) =  0.0;
      t.Sum()(1,0) =  0.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 4:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  0.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 5:
    	t.Sum()(0,0) =  0.0;
      t.Sum()(1,0) =  1.0;
      t.Sum()(2,0) =  1.0;
      return t;
    case 6:
      t.Mult()(0,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (2,0) = -1.0;
      return t;
    case 7:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      t.Sum() (2,0) = -1.0;
      return t;
    case 8:
      t.Mult()(1,0) = -0.5;
      t.Sum() (1,0) =  0.5;
      t.Sum() (2,0) = -1.0;
      return t;
    case 9:
      t.Mult()(2,0) =  1.0;
      return t;
    case 10:
      t.Mult()(2,0) =  1.0;
      t.Sum() (0,0) =  1.0;
      return t;
    case 11:
      t.Mult()(2,0) =  1.0;
      t.Sum() (1,0) =  1.0;
      return t;
    case 12:
      t.Mult()(0,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (2,0) =  1.0;
      return t;
    case 13:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      t.Sum() (2,0) =  1.0;
      return t;
    case 14:
      t.Mult()(1,0) = -0.5;
      t.Sum() (1,0) =  0.5;
      t.Sum() (2,0) =  1.0;
      return t;
    case 15:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Sum() (2,0) = -1.0;
      return t;
    case 16:
      t.Mult()(0,0) =  0.5;
      t.Mult()(2,1) =  1.0;
      t.Sum() (0,0) =  0.5;
      return t;
    case 17:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Mult()(2,1) =  1.0;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 18:
      t.Mult()(1,0) =  0.5;
      t.Mult()(2,1) =  1.0;
      t.Sum() (1,0) =  0.5;
      return t;
    case 19:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Sum() (2,0) =  1.0;
      return t;
    case 20:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      t.Mult()(2,2) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}


int TPZShapePrism::SideNodeLocId(int side, int node)
{
	if(side<6 && node == 0) return side;
	if(side>= 6 && side < 15 && node<2) return SideNodes[side-6][node];
	if(side==15)
	  if (node < 3) return FaceNodes[side-15][node];
	  else if (node == 3) return -1; //previsto para faces triangulares
	
	if(side>15 && side <19 && node <4) return FaceNodes[side-15][node];
	if(side==19)
	  if (node<3) return FaceNodes[side-15][node];
	  else if (node == 3) return -1; // Previsto p/ faces triangulares
	
	if(side==20 && node<6) return node;
	PZError << "TPZShapePrism::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

static int nsidenodes[21] = {
1,1,1,1,1,1,
2,2,2,2,2,2,2,2,2,
3,4,4,4,3,
6};

int TPZShapePrism::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapePrism::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapePrism::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapePrism::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapePrism::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZShapePrism::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapePrism::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapePrism::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

int TPZShapePrism::NSideConnects(int side) {
  if(side<0)   return -1;
  if(side<6)   return 1;//cantos : 0 a 5
  if(side<15)  return 3;//arestas
  if(side==15 || side==19)  return 7;//faces : 15,19 , triangulares
  if(side<19) return 9;//faces : 16 a 18  quadrilaterais
  if(side==20) return 21;//centro : 20
  return -1;
}
//até aqui até aqui até aqui até aqui até aqui até aqui até aqui até aqui até aqui
int TPZShapePrism::SideConnectLocId(int side, int node) {
  if(side<0 || side>20 || node < 0) return -1;
  if(side<6) {
    if(node==0) return side;
  } else
    if(side<9) {//6,7,8
      int s = side-6;//0,1,2
      if(!node) return s;//0,1,2
      if(node==1) return (s+1)%3;//1,2,0
      if(node==2) return side;//6,7,8
    } else
      if(side<12) {//9,10,11
	int s = side-9;//0,1,2
   	if(!node) return s;//0,1,2,4
	if(node==1) return s+3;//3,4,5
	if(node==2) return side;//5,6,7
      } else
	if(side<15) {//12,13,14
	  int s = side-9;//3,4,5
	  if(!node) return s;//3,4,5
	  if(node==1) return (s+1)%3+3;//4,5,3
	  if(node==2) return side;//12,13,14
	} else
	  if(side==15 || side==19) {
	    int s = side-15;
	    if(side==15 && node<7) return FaceConnectLocId[s][node];
	    if(side==19 && node<7) return FaceConnectLocId[s][node];
	  } else
	    if(side<20) {//16,17,18
	      int s = side-15;
	      if(node<9) return FaceConnectLocId[s][node];
	    } else
	      if(side==20 && node<21){
		return node;
	      }
  PZError << "TPZShapePrism::SideConnectLocId called for node = "
	  << node << " and side = " << side << "\n";
  return -1;
}

void TPZShapePrism::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}
