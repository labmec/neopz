/**
 * @file
 * @brief Contains the implementation of the TPZShapeQuad methods.
 */
// $Id: pzshapequad.cpp,v 1.16 2011-05-11 01:47:45 phil Exp $
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
//#include "pzelgq2d.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

using namespace std;

namespace pzshape {
	
	/**Transformation of the point within a quadrilateral face */
	REAL TPZShapeQuad::gTrans2dQ[8][2][2] = {//s* , t*
		{ { 1., 0.},{ 0., 1.} },
		{ { 0., 1.},{ 1., 0.} },
		{ { 0., 1.},{-1., 0.} },
		{ {-1., 0.},{ 0., 1.} },
		{ {-1., 0.},{ 0.,-1.} },//s* = -s   t* = -t  , etc
		{ { 0.,-1.},{-1., 0.} },
		{ { 0.,-1.},{ 1., 0.} },
		{ { 1., 0.},{ 0.,-1.} }
	};
	
	REAL TPZShapeQuad::gRibTrans2dQ1d[4][2] = { {1.,0.},{0.,1.},{-1.,0.},{0.,-1.} };
	
	void TPZShapeQuad::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
		
		REAL x[2],dx[2],y[2],dy[2];
		x[0]  =  (1.-pt[0])/2.;
		x[1]  =  (1.+pt[0])/2.;
		dx[0] = -0.5;
		dx[1] =  0.5;
		y[0]  =  (1.-pt[1])/2.;
		y[1]  =  (1.+pt[1])/2.;
		dy[0] = -0.5;
		dy[1] =  0.5;
		phi(0,0)  = x[0]*y[0];
		phi(1,0)  = x[1]*y[0];
		phi(2,0)  = x[1]*y[1];
		phi(3,0)  = x[0]*y[1];
		dphi(0,0) = dx[0]*y[0];
		dphi(1,0) = x[0]*dy[0];
		dphi(0,1) = dx[1]*y[0];
		dphi(1,1) = x[1]*dy[0];
		dphi(0,2) = dx[1]*y[1];
		dphi(1,2) = x[1]*dy[1];
		dphi(0,3) = dx[0]*y[1];
		dphi(1,3) = x[0]*dy[1];
	}
	
	/*
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input) value of the (4) shape functions
	 * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
	 */
	void TPZShapeQuad::ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi)
	{
		int is;
		for(is=4; is<8; is++)
		{
			phi(is,0) = phi(is%4,0)*phi((is+1)%4,0);
			dphi(0,is) = dphi(0,is%4)*phi((is+1)%4,0)+phi(is%4,0)*dphi(0,(is+1)%4);
			dphi(1,is) = dphi(1,is%4)*phi((is+1)%4,0)+phi(is%4,0)*dphi(1,(is+1)%4);
		}
		phi(8,0) = phi(0,0)*phi(2,0);
		dphi(0,8) = dphi(0,0)*phi(2,0)+phi(0,0)*dphi(0,2);
		dphi(1,8) = dphi(1,0)*phi(2,0)+phi(0,0)*dphi(1,2);

		// Make the generating shape functions linear and unitary
		for(is=4; is<8; is++)
		{
			phi(is,0) += phi(8,0);
			dphi(0,is) += dphi(0,8);
			dphi(1,is) += dphi(1,8);
			phi(is,0) *= 4.;
			dphi(0,is) *= 4.;
			dphi(1,is) *= 4.;
		}
		phi(8,0) *= 16.;
		dphi(0,8) *= 16.;
		dphi(1,8) *= 16.;
	}
	
	void TPZShapeQuad::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
							 TPZFMatrix &phi,TPZFMatrix &dphi) {
		ShapeCorner(pt,phi,dphi);
		int is,d;
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
		REAL out;
		int shape = 4;
		for (int rib = 0; rib < 4; rib++) {
			
			ProjectPoint2dQuadToRib(rib,pt,out);
			TPZVec<int> ids(2);
			TPZManVector<REAL,1> outvec(1,out);
			ids[0] = id[rib%4];
			ids[1] = id[(rib+1)%4];
			REAL store1[20],store2[40];
			int ord2 = order[rib]-1;//two orders : order in x and order in y
			TPZFMatrix phin(ord2,1,store1,20),dphin(2,ord2,store2,40);
			TPZShapeLinear::ShapeInternal(outvec,order[rib],phin,dphin,TPZShapeLinear::GetTransformId1d(ids));
			TransformDerivativeFromRibToQuad(rib,ord2,dphin);
			for (int i = 0; i < ord2; i++) {
				phi(shape,0) = phiblend(rib+4,0)*phin(i,0);
				for(int xj=0;xj<2;xj++) {
					dphi(xj,shape) = dphiblend(xj,rib+4)*phin(i,0)+
					phiblend(rib+4,0)*dphin(xj,i);
				}
				shape++;
			}
		}
		REAL store1[20],store2[40];
		int ord = (order[4]-1)*(order[4]-1);
		TPZFMatrix phin(ord,1,store1,20),dphin(2,ord,store2,40);
		ShapeInternal(pt,order[4]-2,phin,dphin,GetTransformId2dQ(id));
		for(int i=0;i<ord;i++)	{//funcoes de interior s�o em numero ordem-1
			phi(shape,0) = phiblend(8,0)*phin(i,0);
			for(int xj=0;xj<2;xj++) {//x e y
				dphi(xj,shape) = dphiblend(xj,8)*phin(i,0) +
				phiblend(8,0)*dphin(xj,i);
			}
			shape++;
		}
	}
	
	void TPZShapeQuad::SideShape(int side,TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
								 TPZFMatrix &phi,TPZFMatrix &dphi) {
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
				TPZShapePoint::Shape(pt, id, order, phi, dphi);
				break;
			case 4:
			case 5:
			case 6:
			case 7:
				TPZShapeLinear::Shape(pt, id, order, phi, dphi);
				break;
			case 8:
				Shape(pt, id, order, phi, dphi);
		}
	}
	
	void TPZShapeQuad::ShapeInternal(TPZVec<REAL> &x, int order,
									 TPZFMatrix &phi, TPZFMatrix &dphi,int quad_transformation_index) {
		
		if(order < 0) return;
		int ord1 = order+1;
		int numshape = ord1*ord1;
		TPZManVector<REAL> out(2);
		TransformPoint2dQ(quad_transformation_index,x,out);
		
		if(numshape > phi.Rows() || phi.Cols() < 1) phi.Resize(numshape,1);
		if(dphi.Rows() < 2 || dphi.Cols() < numshape) dphi.Resize(2,numshape);
		REAL store1[20],store2[20],store3[20],store4[20];
		TPZFMatrix phi0(ord1,1,store1,20),phi1(ord1,1,store2,20),dphi0(1,ord1,store3,20),dphi1(1,ord1,store4,20);
		
		TPZShapeLinear::fOrthogonal(out[0],ord1,phi0,dphi0);
		TPZShapeLinear::fOrthogonal(out[1],ord1,phi1,dphi1);
		for (int i=0;i<ord1;i++) {
			for (int j=0;j<ord1;j++) {
				int index = i*ord1+j;
				phi(index,0) =  phi0(i,0)* phi1(j,0);
				dphi(0,index) = dphi0(0,i)* phi1(j,0);
				dphi(1,index) =  phi0(i,0)*dphi1(0,j);
			}
		}
		TransformDerivative2dQ(quad_transformation_index,numshape,dphi);
	}
	
	void TPZShapeQuad::TransformDerivative2dQ(int transid, int num, TPZFMatrix &in) {
		
		for(int i=0;i<num;i++) {
			REAL aux[2];
			aux[0] = in(0,i);
			aux[1] = in(1,i);
			in(0,i) = gTrans2dQ[transid][0][0]*aux[0]+gTrans2dQ[transid][1][0]*aux[1];
			in(1,i) = gTrans2dQ[transid][0][1]*aux[0]+gTrans2dQ[transid][1][1]*aux[1];
		}
	}
	
	//transf. o ponto dentro da face quadrilateral
	void TPZShapeQuad::TransformPoint2dQ(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out) {
		
		out[0] = gTrans2dQ[transid][0][0]*in[0]+gTrans2dQ[transid][0][1]*in[1];//Cedric 23/02/99
		out[1] = gTrans2dQ[transid][1][0]*in[0]+gTrans2dQ[transid][1][1]*in[1];//Cedric 23/02/99
	}
	
	void TPZShapeQuad::ProjectPoint2dQuadToRib(int rib, TPZVec<REAL> &in, REAL &out) {
		
		out = gRibTrans2dQ1d[rib][0]*in[0]+gRibTrans2dQ1d[rib][1]*in[1];
	}
	
	int TPZShapeQuad::GetTransformId2dQ(TPZVec<int> &id) {
		
		int id0,id1,minid;
		id0 = (id[0] < id[1]) ? 0 : 1;
		id1 = (id[2] < id[3]) ? 2 : 3;
		minid = (id[id0] < id[id1]) ? id0 : id1;//minid : menor id local
		id0 = (minid+1)%4;//id anterior local (sentido antihorario)
		id1 = (minid+3)%4;//id posterior local (sentido horario)
		minid = id[minid];//minid : menor id global
		
		if (id[id0] < id[id1]) {//antihorario
			
			if (minid == id[0]) return 0;
			if (minid == id[1]) return 2;
			if (minid == id[2]) return 4;
			if (minid == id[3]) return 6;
			
		} else {//horario
			
			if (minid == id[0]) return 1;
			if (minid == id[1]) return 3;
			if (minid == id[2]) return 5;
			if (minid == id[3]) return 7;
		}
		return 0;
	}
	
	void TPZShapeQuad::TransformDerivativeFromRibToQuad(int rib,int num,TPZFMatrix &dphi) {
		
		for (int j = 0;j<num;j++) {
			dphi(1,j) = gRibTrans2dQ1d[rib][1]*dphi(0,j);
			dphi(0,j) = gRibTrans2dQ1d[rib][0]*dphi(0,j);
		}
	}
	
	int TPZShapeQuad::NConnectShapeF(int side, int order) {
		if(side<4) return 1;//0 a 4
		//   int s = side-4;//s = 0 a 14 ou side = 6 a 20
		if(side<8) return (order-1);//6 a 14
		if(side==8) {
			return ((order-1)*(order-1));
		}
		PZError << "TPZShapeQuad::NConnectShapeF, bad parameter side " << side << endl;
		return 0;
	}
	
	int TPZShapeQuad::NShapeF(TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	
	
	
#ifdef _AUTODIFF
	
	void TPZShapeQuad::Shape2dQuadInternal(TPZVec<FADREAL> &x, int order,
										   TPZVec<FADREAL> &phi,int quad_transformation_index) {
		
		const int ndim = 3;
		if(order < 0) return;
		int ord1 = order+1;
		int numshape = ord1*ord1;
		TPZVec<FADREAL> out(2);
		TransformPoint2dQ(quad_transformation_index,x,out);
		
		if(numshape > phi.NElements()/*Rows()*/ || phi[0].size()/*Cols()*/ < ndim) phi.Resize(numshape, FADREAL(ndim, 0.0));
		//if(dphi.Rows() < 2 || dphi.Cols() < numshape) dphi.Resize(2,numshape);
		//REAL store1[20],store2[20],store3[20],store4[20];
		//TPZFMatrix phi0(ord1,1,store1,20),phi1(ord1,1,store2,20),dphi0(1,ord1,store3,20),dphi1(1,ord1,store4,20);
		TPZVec<FADREAL> phi0(20, FADREAL(ndim, 0.0)),
		phi1(20, FADREAL(ndim, 0.0));
		
		TPZShapeLinear::FADfOrthogonal(out[0],ord1,phi0);
		TPZShapeLinear::FADfOrthogonal(out[1],ord1,phi1);
		for (int i=0;i<ord1;i++) {
			for (int j=0;j<ord1;j++) {
				int index = i*ord1+j;
				//phi(index,0) =  phi0(i,0)* phi1(j,0);
				phi[index] =  phi0[i] * phi1[j];
				/*dphi(0,index) = dphi0(0,i)* phi1(j,0);
				 dphi(1,index) =  phi0(i,0)*dphi1(0,j);*/
			}
		}
		//  TransformDerivative2dQ(quad_transformation_index,numshape,phi);
	}
	
	void TPZShapeQuad::TransformPoint2dQ(int transid, TPZVec<FADREAL> &in, TPZVec<FADREAL> &out) {
		
		out[0] = gTrans2dQ[transid][0][0]*in[0]+gTrans2dQ[transid][0][1]*in[1];//Cedric 23/02/99
		out[1] = gTrans2dQ[transid][1][0]*in[0]+gTrans2dQ[transid][1][1]*in[1];//Cedric 23/02/99
	}
	
	/*
	 void TPZShapeQuad::TransformDerivative2dQ(int transid, int num, TPZVec<FADREAL> &in) {
	 
	 for(int i=0;i<num;i++) {
	 double aux[2];
	 aux[0] = in[i].d(0);
	 aux[1] = in[i].d(1);
	 in[i].fastAccessDx(0) = gTrans2dQ[transid][0][0]*aux[0]+gTrans2dQ[transid][1][0]*aux[1];
	 in[i].fastAccessDx(1) = gTrans2dQ[transid][0][1]*aux[0]+gTrans2dQ[transid][1][1]*aux[1];
	 }
	 }
	 
	 void TPZShapeQuad::TransformDerivativeFromRibToQuad(int rib,int num,TPZVec<FADREAL> &phi) {
	 
	 for (int j = 0;j<num;j++) {
	 //dphi(1,j) = gRibTrans2dQ1d[rib][1]*dphi(0,j);
	 //dphi(0,j) = gRibTrans2dQ1d[rib][0]*dphi(0,j);
	 phi[j].fastAccessDx(1) = gRibTrans2dQ1d[rib][1]*phi[j].d(0);
	 phi[j].fastAccessDx(0) = gRibTrans2dQ1d[rib][0]*phi[j].d(0);
	 }
	 }
	 */
	
#endif
	
};



