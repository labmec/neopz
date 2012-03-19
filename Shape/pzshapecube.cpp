/**
 * @file
 * @brief Contains the implementation of the TPZShapeCube methods.
 */
#include "pzshapecube.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
using namespace std;

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
	
//	static int FaceSides[6][4] = { {8,9,10,11},{8,13,16,12},{9,14,17,13},
//		{10,14,18,15},{11,15,19,12},{16,17,18,19} };

	
	void TPZShapeCube::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
		
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
	
	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input) value of the (8) shape functions
	 * @param dphi (input) value of the derivatives of the (8) shape functions holding the derivatives in a column
	 */
	void TPZShapeCube::ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		int is;
		// contribute the ribs
		for(is=8; is<20; is++)
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
		for(is=20; is<26; is++)
		{
			int is1,is2;
			is1 = ContainedSideLocId(is,0);
			is2 = ContainedSideLocId(is,2);
			phi(is,0) = phi(is1,0)*phi(is2,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
			dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
		}
		// contribution of the volume
		for(is=26; is<27; is++)
		{
			int is1,is2;
			is1 = 0;
			is2 = 6;
			phi(is,0) = phi(is1,0)*phi(is2,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
			dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
		}

		// Make the generating shape functions linear and unitary
		for(is=8; is<27; is++)
		{
			TPZStack<int> highsides;
			HigherDimensionSides(is,highsides);
			int h, nh = highsides.NElements();
			for(h=0; h<nh; h++)
			{
				int hs = highsides[h];
				phi(is,0) += phi(hs,0);
				dphi(0,is) += dphi(0,hs);
				dphi(1,is) += dphi(1,hs);
				dphi(2,is) += dphi(2,hs);
			}
			int dim = SideDimension(is);
			int mult = (dim == 1) ? 4 : (dim == 2) ? 16 : (dim ==3) ? 64 : 0;
			phi(is,0) *= mult;
			dphi(0,is) *= mult;
			dphi(1,is) *= mult;
			dphi(2,is) *= mult;
		}

	}
	
	
	void TPZShapeCube::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		ShapeCorner(pt,phi,dphi);
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
			TPZFMatrix<REAL> phin(ordin,1,store1,20),dphin(3,ordin,store2,60);
			phin.Zero();
			dphin.Zero();
			TPZShapeLinear::ShapeInternal(outvalvec,order[rib],phin,dphin,TPZShapeLinear::GetTransformId1d(ids));//ordin = ordem de um lado
			TransformDerivativeFromRibToCube(rib,ordin,dphin);
			for (int i = 0; i < ordin; i++) {
				phi(shape,0) = phiblend(rib+8,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj ,rib+8) *  phin( i, 0) +
					phiblend(rib+8, 0 ) * dphin(xj,i);
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
			TPZFMatrix<REAL> phin(ord,1,store1,20),dphin(3,ord,store2,60);//ponto na face
			phin.Zero();
			dphin.Zero();
			int ordin =  (ord1 > ord2) ? ord1 : ord2;
			ordin--;
			TPZManVector<int> ids(4);
			//TPZVec<int> ids(4);
			//	int id0,id1;
			int i;
			for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
			//id0 = ShapeFaceId[face][0];//numero das shapes da face que compoem a shape atual
			//id1 = ShapeFaceId[face][1];
			TPZShapeQuad::ShapeInternal(outval,ord1-2,phin,dphin,TPZShapeQuad::GetTransformId2dQ(ids));//ordin = ordem de um lado
			TransformDerivativeFromFaceToCube(face,ord,dphin);//ord = numero de shapes
			for(i=0;i<ord;i++)	{
				phi(shape,0) = phiblend(face+20,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj,face+20) * phin(i ,0) +
					phiblend(face+20, 0)* dphin(xj,i);
				}
				shape++;
			}
		}
		//volume shapes
		REAL store1[20],store2[60];
		int ordmin1 = (order[18]-1);
		int ord =  ordmin1*ordmin1*ordmin1;//(p-1)^3 : 0<=n1,n2,n3<=p-2
		TPZFMatrix<REAL> phin(ord,1,store1,20),dphin(3,ord,store2,60);
		phin.Zero();
		dphin.Zero();
		ShapeInternal(pt,ordmin1,phin,dphin);
		for(int i=0;i<ord;i++)	{
			phi(shape,0) = phiblend(26,0)*phin(i,0);
			for(int xj=0;xj<3;xj++) {
				dphi(xj,shape) = dphiblend(xj,26) * phin(i ,0) +
				phiblend(26, 0)*dphin(xj,i);
			}
			shape++;
		}
	}
	
	void TPZShapeCube::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		if(side<0 || side>26) PZError << "TPZCompElC3d::SideShapeFunction. Bad paramenter side.\n";
		else if(side==26) Shape(pt,id,order,phi,dphi);
		else if(side<8) TPZShapePoint::Shape(pt,id,order,phi,dphi);
		else if(side<20) {//8 a 19
			TPZShapeLinear::Shape(pt,id,order,phi,dphi);
		}
		else if(side<26) {//faces
			TPZShapeQuad::Shape(pt,id,order,phi,dphi);
		}
		else
		{
			Shape(pt,id,order,phi,dphi);
		}
	}
	
	void TPZShapeCube::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
									 TPZFMatrix<REAL> &dphi) {//,int cube_transformation_index
		if(order < 1) return;
		int ord = order;//fSideOrder[18]-1;
		order = order*order*order;
		phi.Resize(order,1);
		dphi.Resize(3,order);
		REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
		TPZFMatrix<REAL> phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
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
	
	void TPZShapeCube::TransformDerivativeFromRibToCube(int rib,int num,TPZFMatrix<REAL> &dphi) {
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
	
	void TPZShapeCube::TransformDerivativeFromFaceToCube(int rib,int num,TPZFMatrix<REAL> &dphi) {
		
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
		int in,res = NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	
	/*
	 void TPZShapeCube::LowerDimensionSides(int side,TPZStack<int> &smallsides) {
	 
	 cout << "TPZShapeCube::LowerDimensionSides Nao deve ser usado";
	 DebugStop();
	 
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
	 */
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
			//TPZFMatrix<REAL> phin(ordin,1,store1,20),dphin(3,ordin,store2,60);
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
			//TPZFMatrix<REAL> phin(ord,1,store1,20),dphin(3,ord,store2,60);//ponto na face
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
		//TPZFMatrix<REAL> phin(ord,1,store1,20),dphin(3,ord,store2,60);
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
		//TPZFMatrix<REAL> phi0(ord,1,store1,20),phi1(ord,1,store2,20),phi2(ord,1,store3,20),
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
