// $Id: pzshapeprism.cpp,v 1.7 2007-04-20 18:30:23 caju Exp $
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
using namespace std;

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {

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
  ord = NConnectShapeF(20,order[20-NCornerNodes]);
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

int TPZShapePrism::NShapeF(TPZVec<int> &order) {
  int in,res=NCornerNodes;
  for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
  return res;
}



};
