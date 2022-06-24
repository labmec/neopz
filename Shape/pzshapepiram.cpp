/**
 * @file
 * @brief Contains the implementation of the TPZShapePiram methods.
 */

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
		{ { 1., .5,-.5},{ 0.,-1.,1.}  },//3 2 4 ; original-> 2 3 4 : {-1., .5,-.5},{ 0.,-1.,1.}
		{ {-.5, 1.,-.5},{1.,0.,1.}  } //0 3 4 ; original-> 3 0 4 : {-.5,-1.,-.5},{ 1., 0.,1.}
	};
	
	REAL TPZShapePiram::gFaceSum3dPiram2d[5][2] = { {.0,0.},{-.5,0.},{-.5,0.},{-.5,0.},{-.5,0.} };//{ {.0,0.},{-.5,0.},{-.5,0.},{-.5,0.},{-.5,0.} };//original
	
	
	
    /**
     * Computes the generating shape functions for a quadrilateral element
     * @param pt (input) point where the shape function is computed
     * @param phi (input/output) value of the (4) shape functions
     * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
     */
    void TPZShapePiram::ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
    {
        int is;
        // contribute the ribs
        for(is=NCornerNodes; is<NCornerNodes+8; is++)
        {
            if(nshape[is-NCornerNodes] < 1) continue;
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
            if(nshape[is-NCornerNodes] < 1) continue;
            int is1,is2,is3;
            is1 = ContainedSideLocId(is,0);
            is2 = ContainedSideLocId(is,1);
            is3 = ContainedSideLocId(is,2);
            phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
            dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
            dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
            dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
            
        }
        if(nshape[NSides-NCornerNodes-1] > 0)
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
                if(nshape[is-NCornerNodes] < 1) continue;
                phi(is,0) *= sidescale[is];
                dphi(0,is) *= sidescale[is];
                dphi(1,is) *= sidescale[is];
                dphi(2,is) *= sidescale[is];
            }
        }
        
        
    }
    
	void TPZShapePiram::Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
		CornerShape(pt,phi,dphi);
		bool linear = true;
		int is,d;
        for(is=NCornerNodes; is<NSides; is++){
            if(order[is-NCornerNodes] > 1){
                linear = false;
            }
        }
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
			TPZVec<int64_t> ids(2);
			TPZManVector<REAL,1> outvalvec(1,outval);
			int id0,id1;
			id0 = SideNodes[rib][0];
			id1 = SideNodes[rib][1];
			ids[0] = id[id0];
			ids[1] = id[id1];
			REAL store1[20],store2[60];
			int nshape = order[rib]-1;//three orders : order in x , order in y and order in z
			TPZFMatrix<REAL> phin(nshape,1,store1,20),dphin(3,nshape,store2,60);
			phin.Zero();
			dphin.Zero();
			TPZShapeLinear::ShapeInternal(outvalvec,order[rib],phin,dphin,TPZShapeLinear::GetTransformId1d(ids));//ordin = ordem de um lado
			TransformDerivativeFromRibToPiram(rib,nshape,dphin);
			for (int i = 0; i < nshape; i++) {
				phi(shape,0) = phiblend(rib+5,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj ,rib+5) * phin( i, 0) +
					phiblend(rib+5, 0 )  * dphin(xj,i);
				}
				shape++;
			}
		}

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
			int nshape;
			if(!face) nshape = (ord1-1)*(ord1-1);//face quadrada
			else nshape = (ord1-2)*(ord1-1)/2;//face triangular
			TPZFNMatrix<60> phin(nshape,1),dphin(3,nshape);//ponto na face
			phin.Zero();
			dphin.Zero();
			TPZManVector<int64_t> ids(4);
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
			TransformDerivativeFromFaceToPiram(face,nshape,dphin);//ordin = numero de shapes
			for(i=0;i<nshape;i++)	{
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
		int nshape=0,i;
		for(i=0;i<order[13]-2;i++) {
			nshape += (i+1)*(i+2) / 2;
		}
		TPZFNMatrix<40> phin(nshape,1),dphin(3,nshape);
		phin.Zero();
		dphin.Zero();
		ShapeInternal(pt,order[13],phin,dphin);
		for(i=0;i<nshape;i++)	{
			phi(shape,0) = phiblend(NSides-1,0)*phin(i,0);
			for(int xj=0;xj<3;xj++) {
				dphi(xj,shape) = dphiblend(xj,NSides-1) * phin(i ,0) +
				phiblend(NSides-1, 0) * dphin(xj,i);
			}
			shape++;
		}
	}
	
	void TPZShapePiram::SideShape(int side, TPZVec<REAL> &point, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
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
    
    void TPZShapePiram::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        //DebugStop();
        int64_t nsides = TPZShapePiram::NSides;
        int nshape;
        
        int linha = 0;
        for (int side = 0; side < nsides; side++)
        {
            
            nshape = 1;
            if(side >= NCornerNodes) nshape = NConnectShapeF(side,order[side-NCornerNodes]);
            int sideorder = 1;
            if(side >= NCornerNodes) sideorder = order[side-NCornerNodes];
            
            TPZGenMatrix<int> locshapeorders(nshape,3);
            SideShapeOrder(side, id, sideorder, locshapeorders);
            
            int nlin = locshapeorders.Rows();
            int ncol = locshapeorders.Cols();
            
            for (int il = 0; il<nlin; il++)
            {
                for (int jc = 0; jc<ncol; jc++)
                {
                    shapeorders(linha, jc) = locshapeorders(il, jc);
                }
                linha++;
            }
        }
    }
    
    
    void TPZShapePiram::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        //DebugStop();
        if (side<=4)
        {
            if (shapeorders.Rows() != 1)
            {
                DebugStop();
            }
            shapeorders(0,0) = 1;
            shapeorders(0,1) = 0;
            shapeorders(0,2) = 0;
        }
        else if (side>=5 && side<=12)
        {
            int nshape = order-1;
            if (shapeorders.Rows() != nshape)
            {
                DebugStop();
            }
            for (int ioy = 0; ioy < order-1; ioy++)
            {
                shapeorders(ioy,0) = ioy+2;
            }
        }
        else if (side == 13)
        {
            // quadrilatero
            if (shapeorders.Rows() != (order-1)*(order-1))
            {
                DebugStop();
            }
            TPZStack<int> lowersides;
            LowerDimensionSides(side, lowersides);
            lowersides.Push(side);
            
            int nnodes = NSideNodes(side);
            
            TPZManVector<int64_t, 4> locid(nnodes);
            for (int node=0; node<locid.size(); node++) {
                locid[node] = id[ContainedSideLocId(side, node)];// SideNodeLocId( side, node);
            }
            
            int nshape = (order-1)*(order-1);
            TPZGenMatrix<int> locshapeorders(nshape,3);
            
            
            TPZShapeQuad::SideShapeOrder(8,locid, order, locshapeorders);
            
            // temos que arrumar a saida de locshapeorders para adequar a orientacao dos vetores que geram
            // a face do lado side
            
            // aqui o locshapeorders esta armazenado so para x e y
            for (int il = 0; il<nshape; il++)
            {
                shapeorders(il, 0) = locshapeorders(il, 0);
                shapeorders(il, 1) = locshapeorders(il, 1);
                shapeorders(il, 2) = locshapeorders(il, 2);
            }
        }
        else if (side >= 14 && side <= 17)
        {
            //triangulos
            int nshape = (order-2)*(order-1)/2;
            if (shapeorders.Rows() != nshape)
            {
                DebugStop();
            }
            TPZStack<int> lowersides;
            LowerDimensionSides(side, lowersides);
            lowersides.Push(side);
            
            int nnodes = NSideNodes(side);
            
            TPZManVector<int64_t, 4> locid(nnodes);
            for (int node=0; node<locid.size(); node++) {
                locid[node] = id[ContainedSideLocId(side, node)];
            }
            
            TPZGenMatrix<int> locshapeorders(nshape,3);
            
            
            TPZShapeTriang::SideShapeOrder(6,locid, order, locshapeorders);
            
            for (int il = 0; il<nshape; il++)
            {
                shapeorders(il, 0) = locshapeorders(il, 0);
                shapeorders(il, 1) = locshapeorders(il, 1);
                shapeorders(il, 2) = locshapeorders(il, 2);
            }
        }
        else
        {
            // interno
            int totsum = 0;
            for(int i=0;i<order - 2;i++) {
                totsum += (i+1)*(i+2) / 2;
            }
            int nshape = totsum;
            
            if (shapeorders.Rows() != nshape) {
                DebugStop();
            }
            int count = 0;
            int ord = order-2;
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        int soma = i+j+k;
                        if( i+j+k < ord ) // Duvida
                        {
                            shapeorders(count,0) = 3 + soma;
                            shapeorders(count,1) = 3 + soma;
                            shapeorders(count,2) = 3 + soma;
                            count++;
                        }
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
	
	int TPZShapePiram::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}


    int TPZShapePiram::ClassId() const{
        return Hash("TPZShapePiram") ^ pztopology::TPZPyramid::ClassId() << 1;
    }
	
};
