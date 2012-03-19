/**
 * @file
 * @brief Contains the implementation of the TPZShapePiram methods.
 */
// $Id: pzshapepiram.cpp,v 1.12 2011-05-11 01:47:45 phil Exp $
#include "pzshapepiram.h"
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
	REAL TPZShapePiram::gRibTrans3dPiram1d[8][3] = {//par�metros de arestas
		{ 1., 0.,0.} , { 0., 1.,0.} ,//percorre o sentido
		{-1., 0.,0.} , { 0.,-1.,0.} ,//da aresta segundo
		{ .5, .5,1.} , {-.5,.5 ,1.} ,//SideNodes[8][2]
		{-.5,-.5,1.} , { .5,-.5,1.}
	};
	//n�o tem vetor associado -> OK!
	
	//Projection of the point within a piramide to a face
	REAL TPZShapePiram::gFaceTrans3dPiram2d[5][2][3] = {//par�metros de faces
		{ { 1., 0., 0.},{ 0., 1.,0.}  },//0 1 2 3  : percorre os eixos segundo
		{ { 1.,-.5,-.5},{ 0., 1.,1.}  },//0 1 4    : FaceNodes[5][4]
		{ { .5, 1.,-.5},{-1., 0.,1.}  },//1 2 4
		//{ {-1., .5,-.5},{ 0.,-1.,1.}  },//ORIGINAL
		//{ {-.5,-1.,-.5},{ 1., 0.,1.}  }//ORIGINAL
		{ { 1., .5,-.5},{ 0.,-1.,1.}  },//3 2 4 ; original-> 2 3 4 : {-1., .5,-.5},{ 0.,-1.,1.}
		{ {-.5, 1.,-.5},{1.,0.,1.}  } //0 3 4 ; original-> 3 0 4 : {-.5,-1.,-.5},{ 1., 0.,1.}
	};
	
	REAL TPZShapePiram::gFaceSum3dPiram2d[5][2] = { {.0,0.},{-.5,0.},{-.5,0.},{-.5,0.},{-.5,0.} };//{ {.0,0.},{-.5,0.},{-.5,0.},{-.5,0.},{-.5,0.} };//original
	
	void TPZShapePiram::CornerShape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
		
		/*if(abs(pt[0])<1.e-10 && abs(pt[1])<1.e-10 && pt[2]==1.) {
		 //para testes com transforma��es geometricas
		 //(0,0,1) nunca � um ponto de integra��o
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
		REAL T0xz,T0yz,T1xz,T1yz;
		if (fabs(1.-pt[2]) < 1.e-8) 
		{
			if (fabs(pt[0]) > 1.e-8 || fabs(pt[1]) > 1.e-8) 
			{
				DebugStop();
			}
			T0xz = 0.5;
			T0yz = 0.5;
			T1xz = 0.5;
			T1yz = 0.5;
		}
		else 
		{
			T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
			T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
			T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
			T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
		}
		
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
	
	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input/output) value of the (4) shape functions
	 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
	 */
	void TPZShapePiram::ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		int is;
		// contribute the ribs
		for(is=NCornerNodes; is<NCornerNodes+8; is++)
		{
			int is1,is2;
			is1 = ContainedSideLocId(is,0);
			is2 = ContainedSideLocId(is,1);
			phi(is,0) = phi(is1,0)*phi(is2,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
			dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
		}
		// contribution of the faces
		// quadrilateral face
		is = 13;
		{
			int is1,is2;
			is1 = ShapeFaceId[0][0];// SideConnectLocId(is,0);
			is2 = ShapeFaceId[0][2];// SideConnectLocId(is,2);
			phi(is,0) = phi(is1,0)*phi(is2,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
			dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
		}
		is++;
		// triangular faces
		for(;is<18; is++)
		{
			int is1,is2,is3;
			is1 = ContainedSideLocId(is,0);
			is2 = ContainedSideLocId(is,1);
			is3 = ContainedSideLocId(is,2);
			phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
			dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
			
		}
		{
			int is1 = 0;
			int is2 = 2;
			int is3 = 4;
			phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
			dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
		}
		

		// Make the generating shape functions linear and unitary
		// contribute the ribs
		for(is=NCornerNodes; is<NCornerNodes+4; is++)
		{
			int isface = 13;
			phi(is,0) += phi(isface,0);
			dphi(0,is) += dphi(0,isface);
			dphi(1,is) += dphi(1,isface);
			dphi(2,is) += dphi(2,isface);
		}
		// scaling the shapefunctions
		{
			REAL sidescale[] = {1.,1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,4.,4.,16.,27.,27.,27.,27.,64.};
			for(is=5; is<NSides; is++)
			{
				phi(is,0) *= sidescale[is];
				dphi(0,is) *= sidescale[is];
				dphi(1,is) *= sidescale[is];
				dphi(2,is) *= sidescale[is];
			}
		}

		
	}
	
	/*REAL TPZCompEl::T(REAL &pt0,REAL &pt1) {
	 return ( ( ( 1 - pt1 ) - pt0 ) / 2 / (1 - pt1) );
	 }
	 
	 void TPZCompEl::ShapeCornerPira(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	 
	 
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
	
	void TPZShapePiram::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
		CornerShape(pt,phi,dphi);
		bool linear = true;
		int is,d;
		for(is=NCornerNodes; is<NSides; is++) if(order[is-NCornerNodes] > 1) linear = false;
		if(linear) return;
		
		TPZFNMatrix<100> phiblend(NSides,1),dphiblend(Dimension,NSides);
		for(is=0; is<NCornerNodes; is++)
		{
			phiblend(is,0) = phi(is,0);
			for(d=0; d<Dimension; d++)
			{
				dphiblend(d,is) = dphi(d,is);
			}
		}
		ShapeGenerating(pt,phiblend,dphiblend);
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
			TPZFMatrix<REAL> phin(ordin,1,store1,20),dphin(3,ordin,store2,60);
			phin.Zero();
			dphin.Zero();
			TPZShapeLinear::ShapeInternal(outvalvec,order[rib],phin,dphin,TPZShapeLinear::GetTransformId1d(ids));//ordin = ordem de um lado
			TransformDerivativeFromRibToPiram(rib,ordin,dphin);
			for (int i = 0; i < ordin; i++) {
				phi(shape,0) = phiblend(rib+5,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj ,rib+5) * phin( i, 0) +
					phiblend(rib+5, 0 )  * dphin(xj,i);
				}
				shape++;
			}
		}
		//if(order[13]<2) return;//ordem do elemento
		//face shapes
		for (int face = 0; face < 5; face++) {
			if (order[face+8] < 2) continue;
			if(face>0 && order[face+8]<=2) continue;//s� a face 13 tem shape associada com ordem p=2
			TPZManVector<REAL,2> outval(2);
			ProjectPoint3dPiramToFace(face,pt,outval);
			int ord1;//,ord2;
			ord1 = order[8+face];
			//ord2 = ord1;
			//FaceOrder(face,ord1,ord2);//ordem da face
			if(face && ord1<3) continue;//uma face com ordem < 3 n�o tem shape associada
			int ordin;
			if(!face) ordin = (ord1-1)*(ord1-1);//face quadrada
			else ordin = (ord1-2)*(ord1-1)/2;//face triangular
			TPZFNMatrix<60> phin(ordin,1),dphin(3,ordin);//ponto na face
			phin.Zero();
			dphin.Zero();
			TPZManVector<int> ids(4);
			//	int id0,id1,id2;
			int i;
			if(!face) for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
			else for(i=0;i<3;i++) ids[i] = id[FaceNodes[face][i]];
			//	id0 = ShapeFaceId[face][0];//indice das shapes da face que compoem a shape atual
			//	id1 = ShapeFaceId[face][1];//equivale a FaceIdsCube(face,ids,id,id0,id1);
			// 	id2 = ShapeFaceId[face][2];//if(face == 0) id3 = ShapeFaceId[face][3];
			int transid;
			if(!face) {
				transid = TPZShapeQuad::GetTransformId2dQ(ids);
				TPZShapeQuad::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
			} else {
				ids.Resize(3);
				transid = TPZShapeTriang::GetTransformId2dT(ids);
				outval[0] = (outval[0]+1.)/2.;//devido a corre��o na fun��o
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
				phi(shape,0) = phiblend(face+13,0)*phin(i,0);//face quadril�teral
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj,face+13)* phin(i ,0) +
					phiblend(face+13, 0)*dphin(xj,i);
				}
				shape++;
			}
		}
		if(order[13]<3) return;//n�o h� ordens para cantos
		//volume shapes
		int ord=0,i;
		for(i=0;i<order[13]-2;i++) {
			ord += (i+1)*(i+2) / 2;
		}
		TPZFNMatrix<40> phin(ord,1),dphin(3,ord);
		phin.Zero();
		dphin.Zero();
		ShapeInternal(pt,order[13],phin,dphin);
		for(i=0;i<ord;i++)	{
			phi(shape,0) = phiblend(NSides-1,0)*phin(i,0);
			for(int xj=0;xj<3;xj++) {
				dphi(xj,shape) = dphiblend(xj,NSides-1) * phin(i ,0) +
				phiblend(NSides-1, 0) * dphin(xj,i);
			}
			shape++;
		}
	}
	
	void TPZShapePiram::SideShape(int side, TPZVec<REAL> &point, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
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
	
	void TPZShapePiram::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
									  TPZFMatrix<REAL> &dphi) {
		//valor da fun��o e derivada do produto das fun��es ortogonais
		if(order < 3) return;
		int ord = order-2;
		REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
		TPZFMatrix<REAL> phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
		dphi0(1,ord,store4,20),dphi1(1,ord,store5,20),dphi2(1,ord,store6,20);
		TPZShapeLinear::fOrthogonal(x[0],ord,phi0,dphi0);//f e df            -1<=x0<=1
		TPZShapeLinear::fOrthogonal(x[1],ord,phi1,dphi1);//g e dg            -1<=x1<=1
		TPZShapeLinear::fOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);//h e dh       0<=x3<=1 -> -1<=2*x2-1<=1
		int index = 0;//x � ponto de integra��o dentro da pir�mide
		for (int i=0;i<ord;i++) {
			for (int j=0;j<ord;j++) {
				for (int k=0;k<ord;k++) {
					if( i+j+k < ord ) {
						//int index = ord*(ord*i+j)+k; //i,j,k � o grau das fun��es ortogonais
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
	
	void TPZShapePiram::TransformDerivativeFromRibToPiram(int rib,int num,TPZFMatrix<REAL> &dphi) {
		for (int j = 0;j<num;j++) {
			dphi(2,j) = gRibTrans3dPiram1d[rib][2]*dphi(0,j);
			dphi(1,j) = gRibTrans3dPiram1d[rib][1]*dphi(0,j);
			dphi(0,j) = gRibTrans3dPiram1d[rib][0]*dphi(0,j);
		}
	}
	
	void TPZShapePiram::TransformDerivativeFromFaceToPiram(int face,int num,TPZFMatrix<REAL> &dphi) {
		for (int j = 0;j<num;j++) {
			dphi(2,j) = gFaceTrans3dPiram2d[face][0][2]*dphi(0,j)+gFaceTrans3dPiram2d[face][1][2]*dphi(1,j);
			REAL dphi1j = dphi(1,j);
			dphi(1,j) = gFaceTrans3dPiram2d[face][0][1]*dphi(0,j)+gFaceTrans3dPiram2d[face][1][1]*dphi(1,j);
			dphi(0,j) = gFaceTrans3dPiram2d[face][0][0]*dphi(0,j)+gFaceTrans3dPiram2d[face][1][0]*dphi1j;//dphi(1,j);
		}
	}
	//transforma a derivada no ponto dentro da face
	void TPZShapePiram::TransformDerivativeFace3dPiram(int transid, int face, int num, TPZFMatrix<REAL> &in) {
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
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	
	
};
