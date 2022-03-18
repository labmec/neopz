/**
 * @file
 * @brief Contains the implementation of the TPZShapePrism methods.
 */

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
	REAL TPZShapePrism::gRibTrans3dPrisma1d[9][3] = {//par�metros de arestas
		{ 2., 1.,0.} , {-1., 1.,0.} ,//percorre o sentido
		{-1.,-2.,0.} , { 0., 0.,1.} ,//da aresta segundo : { 1., 2.,0.} , { 0., 0.,1.}
		{ 0., 0.,1.} , { 0., 0.,1.} ,//SideNodes[9][2]
		{ 2., 1.,0.} , {-1., 1.,0.} ,
		{-1.,-2.,0.}                                     //{ 1., 2.,0.}
	};
	
	REAL TPZShapePrism::gRibSum3dPrisma1d[9] = {-1.,0.,1.,0.,0.,0.,-1.,0.,1.};//{-1.,0.,-1.,0.,0.,0.,-1.,0.,-1.};
	
	//Projection of the point within a piramide to a face
	REAL TPZShapePrism::gFaceTrans3dPrisma2d[5][2][3] = {//par�metros de faces
		{ { 2., 0., 0.},{ 0., 2., 0.}  },//0 1 2   : percorre os eixos segundo
		{ { 2., 0., 0.},{ 0., 0., 1.}  },//0 1 4 3    : FaceNodes[5][4]
		{ {-1., 1., 0.},{ 0., 0., 1.}  },//1 2 5 4
		{ { 0., 2., 0.},{ 0., 0., 1.}  },//0 2 5 3
		{ { 2., 0., 0.},{ 0., 2., 0.}  } //3 4 5
	};
	
	REAL TPZShapePrism::gFaceSum3dPrisma2d[5][2] = { {-1.,-1.},{-1.,0.},{0.,0.},{-1.,0.},{-1.,-1.} };
	
    void TPZShapePrism::ShapeCorner(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
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
    //esse metodo nao vai ser mais utilizado
	void TPZShapePrism::CornerShape(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
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
	
	
	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input) value of the (4) shape functions
	 * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
	 */
	void TPZShapePrism::ShapeGenerating(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		int is;
		// 9 ribs
		for(is=6; is<NSides; is++)
		{
			int nsnodes = NSideNodes(is);
			switch(nsnodes)
			{
				case 2:
				{
					int is1 = SideNodeLocId(is,0);
					int is2 = SideNodeLocId(is,1);
					phi(is,0) = phi(is1,0)*phi(is2,0);
					dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
					dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
					dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
				}
					break;
				case 3:
				{
					int is0 = 0;
					int is1 = 1;
					int is2 = 2;
					int is3 = 3;
					int is4 = 4;
					if(is == 19)
					{
						is2 = 5;
					}
					phi(is,0) = (phi(is0,0)+phi(is3,0))*(phi(is1,0)+phi(is4,0))*phi(is2,0);
					int d;
					for(d=0; d<3; d++)
					{
						dphi(d,is) = 
						(dphi(d,is0)+dphi(d,is3))*(phi(is1,0)+phi(is4,0))*phi(is2,0) +
						(phi(is0,0)+phi(is3,0))*(dphi(d,is1)+dphi(d,is4))*phi(is2,0) +
						(phi(is0,0)+phi(is3,0))*(phi(is1,0)+phi(is4,0))*dphi(d,is2);
					}
				}
					break;
				case 4:
				{
					int is1 = SideNodeLocId(is,0);
					int is2 = SideNodeLocId(is,2);
					phi(is,0) = phi(is1,0)*phi(is2,0);
					dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
					dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
					dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
				}
					break;
				case 6:
				{
					int is1 = 0;
					int is2 = 4;
					int is3 = 2;
					int is4 = 5;
					phi(is,0) = phi(is1,0)*phi(is2,0)*(phi(is3,0)+phi(is4,0));
					dphi(0,is) = dphi(0,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(0,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(0,is3)+dphi(0,is4));
					dphi(1,is) = dphi(1,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(1,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(1,is3)+dphi(1,is4));
					dphi(2,is) = dphi(2,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(2,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(2,is3)+dphi(2,is4));
				}
					break;
				default:
					;
			}
		}
		// Make the generating shape functions linear and unitary
		for(is=6; is<NSides; is++)
		{
			TPZStack<int> highsides;
			HigherDimensionSides(is,highsides);
			int h, nh = highsides.NElements();
			for(h=0; h<nh; h++)
			{
				int hs = highsides[h];
				if(NSideNodes(hs) != 4) continue;
				phi(is,0) += phi(hs,0);
				dphi(0,is) += dphi(0,hs);
				dphi(1,is) += dphi(1,hs);
				dphi(2,is) += dphi(2,hs);
			}
		}
		REAL mult[] = {1.,1.,1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,4.,4.,4.,27.,16.,16.,16.,27.,8.};
		for(is=6;is<NSides; is++)
		{
			phi(is,0) *= mult[is];
			dphi(0,is) *= mult[is];
			dphi(1,is) *= mult[is];
			dphi(2,is) *= mult[is];
		}
		
	}
	
    /**
     * Computes the generating shape functions for a quadrilateral element
     * @param pt (input) point where the shape function is computed
     * @param phi (input) value of the (4) shape functions
     * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
     */
    void TPZShapePrism::ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
    {
        int is;
        // 9 ribs
        for(is=6; is<NSides; is++)
        {
            int nsnodes = NSideNodes(is);
            switch(nsnodes)
            {
                case 2:
                {
                    if(nshape[is-6] < 1) continue;
                    int is1 = SideNodeLocId(is,0);
                    int is2 = SideNodeLocId(is,1);
                    phi(is,0) = phi(is1,0)*phi(is2,0);
                    dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                    dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                    dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
                }
                    break;
                case 3:
                {
                    if(nshape[is-6] < 1) continue;
                    int is0 = 0;
                    int is1 = 1;
                    int is2 = 2;
                    int is3 = 3;
                    int is4 = 4;
                    if(is == 19)
                    {
                        is2 = 5;
                    }
                    phi(is,0) = (phi(is0,0)+phi(is3,0))*(phi(is1,0)+phi(is4,0))*phi(is2,0);
                    int d;
                    for(d=0; d<3; d++)
                    {
                        dphi(d,is) =
                        (dphi(d,is0)+dphi(d,is3))*(phi(is1,0)+phi(is4,0))*phi(is2,0) +
                        (phi(is0,0)+phi(is3,0))*(dphi(d,is1)+dphi(d,is4))*phi(is2,0) +
                        (phi(is0,0)+phi(is3,0))*(phi(is1,0)+phi(is4,0))*dphi(d,is2);
                    }
                }
                    break;
                case 4:
                {
                    int is1 = SideNodeLocId(is,0);
                    int is2 = SideNodeLocId(is,2);
                    phi(is,0) = phi(is1,0)*phi(is2,0);
                    dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                    dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                    dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
                }
                    break;
                case 6:
                {
                    if(nshape[is-6] < 1) continue;
                    int is1 = 0;
                    int is2 = 4;
                    int is3 = 2;
                    int is4 = 5;
                    phi(is,0) = phi(is1,0)*phi(is2,0)*(phi(is3,0)+phi(is4,0));
                    dphi(0,is) = dphi(0,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(0,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(0,is3)+dphi(0,is4));
                    dphi(1,is) = dphi(1,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(1,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(1,is3)+dphi(1,is4));
                    dphi(2,is) = dphi(2,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(2,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(2,is3)+dphi(2,is4));
                }
                    break;
                default:
                    ;
            }
        }
        // Make the generating shape functions linear and unitary
        for(is=6; is<NSides; is++)
        {
            if(nshape[is-6] < 1) continue;
            TPZStack<int> highsides;
            HigherDimensionSides(is,highsides);
            int h, nh = highsides.NElements();
            for(h=0; h<nh; h++)
            {
                int hs = highsides[h];
                if(NSideNodes(hs) != 4) continue;
                phi(is,0) += phi(hs,0);
                dphi(0,is) += dphi(0,hs);
                dphi(1,is) += dphi(1,hs);
                dphi(2,is) += dphi(2,hs);
            }
        }
        REAL mult[] = {1.,1.,1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,4.,4.,4.,27.,16.,16.,16.,27.,8.};
        for(is=6;is<NSides; is++)
        {
            if(nshape[is-6] < 1) continue;
            phi(is,0) *= mult[is];
            dphi(0,is) *= mult[is];
            dphi(1,is) *= mult[is];
            dphi(2,is) *= mult[is];
        }
        
    }
	void TPZShapePrism::Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
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
		//  if(order[14]<2) return;//order tem as ordens dos lados do elemento
		int shape = 6;
		//rib shapes
		for (int rib = 0; rib < 9; rib++) {//todas as arestas
			REAL outval;
			ProjectPoint3dPrismaToRib(rib,pt,outval);
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
			TransformDerivativeFromRibToPrisma(rib,nshape,dphin);
			for (int i = 0; i < nshape; i++) {
				phi(shape,0) = phiblend(rib+6,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj ,rib+6)  * phin( i, 0) +
					phiblend(rib+6, 0 )   * dphin(xj,i);
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
			if((face==0 || face==4) && ord1<3) continue;//uma face triangular com ordem < 3 n�o tem shape associada
			int ordin;
			if(face && face<4) ordin = (ord1-1)*(ord1-1);//faces quadrilaterais
			else ordin = (ord1-2)*(ord1-1)/2;//face triangular
			TPZFMatrix<REAL> phin(ordin,1,store1,20),dphin(3,ordin,store2,60);//ponto na face
			phin.Zero();
			dphin.Zero();
			TPZManVector<int64_t> ids(4);
			//	int id0,id1,id2
			int i;
			if(face ==0 || face == 4) for(i=0;i<3;i++) ids[i] = id[FaceNodes[face][i]];
			else for(i=0;i<4;i++) ids[i] = id[FaceNodes[face][i]];
			//    id0 = ShapeFaceId[face][0];//indice das shapes da face x
			//    id1 = ShapeFaceId[face][1];//que compoem a shape atual
			//    id2 = ShapeFaceId[face][2];
			int transid;
			if(face && face<4) {
				transid = TPZShapeQuad::GetTransformId2dQ(ids);
				TPZShapeQuad::ShapeInternal(outval,ord1-2,phin,dphin,transid);//ordin = ordem de um lado
			} else {
				ids.Resize(3);
				transid = TPZShapeTriang::GetTransformId2dT(ids);
				outval[0] = (outval[0]+1.)/2.;//devido a corre��o na fun��o
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
				phi(shape,0) = phiblend(face+15,0)*phin(i,0);//face quadril�teral
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj,face+15) * phin(i ,0) +
					phiblend(face+15, 0) * dphin(xj,i);
				}
				shape++;
			}
		}
		if(order[14]<3) return;
		//volume shapes
		int ord=0;
		ord = NConnectShapeF(20,order[20-NCornerNodes]);
		TPZFNMatrix<60> phin(ord,1),dphin(3,ord);
		phin.Zero();
		dphin.Zero();
		ShapeInternal(pt,order[14],phin,dphin);
		for(int i=0;i<ord;i++)	{
			phi(shape,0) = phiblend(NSides-1,0)*phin(i,0);
			for(int xj=0;xj<3;xj++) {
				dphi(xj,shape) = dphiblend(xj,NSides-1)* phin(i ,0) +
				phiblend(NSides-1, 0)* dphin(xj,i);
			}
			shape++;
		}
		
		//}//while
		
	}
	
	void TPZShapePrism::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
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
    
    void TPZShapePrism::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        //DebugStop();
        
        int64_t nsides = TPZShapePrism::NSides;
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
    
    
    void TPZShapePrism::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        //DebugStop();
        if (side<=5)
        {
            if (shapeorders.Rows() != 1)
            {
                DebugStop();
            }
            shapeorders(0,0) = 1;
            shapeorders(0,1) = 0;
            shapeorders(0,2) = 0;
        }
        else if (side>5 && side<15)
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
        else if (side == 15||side == 19)
        {
            int nshape = (order-2)*(order-1)/2;
            if (shapeorders.Rows() != nshape)
            {
                DebugStop();
            }
            TPZStack<int> lowersides;
            LowerDimensionSides(side, lowersides);
            lowersides.Push(side);
            
            //TPZVec<int> locsideorder(lowersides.size(),order);
            
            int nnodes = NSideNodes(side);
            
            TPZManVector<int64_t, 4> locid(nnodes);
            for (int node=0; node<locid.size(); node++) {
                locid[node] = id[ContainedSideLocId(side, node)];// SideNodeLocId( side, node);
            }// sera que esta pegando os ids corretos mesmo?
            
            TPZGenMatrix<int> locshapeorders(nshape,3);
            
            
            TPZShapeTriang::SideShapeOrder(6,locid, order, locshapeorders);
            
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
        else if (side >= 16 && side <=18)
        {
            if (shapeorders.Rows() != (order-1)*(order-1))
            {
                DebugStop();
            }
            TPZStack<int> lowersides;
            LowerDimensionSides(side, lowersides);
            lowersides.Push(side);
            
            //TPZVec<int> locsideorder(lowersides.size(),order);
            
            int nnodes = NSideNodes(side);
            
            TPZManVector<int64_t, 4> locid(nnodes);
            for (int node=0; node<locid.size(); node++) {
                locid[node] = id[ContainedSideLocId(side, node)];// SideNodeLocId( side, node);
            }// sera que esta pegando os ids corretos mesmo?
            
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
        else
        {   // interno
            int nshape = (order-2)*(order-1)*(order-1)/2;
            if (shapeorders.Rows() != nshape) {
                DebugStop();
            }
            int count = 0;
            int ord1 = order - 2;
            int ord2 = order - 1;
            for (int i=0; i<ord1; i++) {
                for (int j=0; j<ord1; j++) {
                    for (int k=0; k<ord2; k++) {
                        int a = i;
                        int b = j;
                        int c = k;
                        int maxAB = a+b;//a>b? a : b;
                        if (   ( (a+b)<ord1 )  && (c < ord2)   ) // Duvida
                        {
                            shapeorders(count,0) = 3 + maxAB;
                            shapeorders(count,1) = 3 + maxAB;
                            shapeorders(count,2) = 2 + c;
                            count++;
                        }
                        
                    }
                    
                }
            }

//            int orderplus1 = order + 2;
//            int orderplus2 = order + 1;
//            for (int i=3; i<orderplus1; i++) {
//                for (int j=3; j<orderplus1; j++) {
//                    for (int k=3; k<orderplus2; k++) {
//                        int a = i;
//                        int b = j;
//                        int c = k;
//                        if (   ( (a+b)<orderplus1 )  && (c < orderplus2)   ) // Duvida
//                        {
//                            shapeorders(count,0) = a;
//                            shapeorders(count,1) = b;
//                            shapeorders(count,2) = c;
//                            count++;
//                        }
//                        
//                    }
//                    
//                }
//            }

        }

    }
	
	
	void TPZShapePrism::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
									  TPZFMatrix<REAL> &dphi) {
		//valor da fun��o e derivada do produto das fun��es ortogonais
		if(order < 3) return;
		int ord1 = order-2;
		int ord2 = order-1;
		
		TPZFNMatrix<20,REAL> phi0(ord1,1),phi1(ord1,1),phi2(ord2,1),
		dphi0(1,ord1),dphi1(1,ord1),dphi2(1,ord2);
		TPZShapeLinear::fOrthogonal(2.*x[0]-1.,ord1,phi0,dphi0);//f e df       0<=x0<=1 -> -1<=2*x0-1<=1
		TPZShapeLinear::fOrthogonal(2.*x[1]-1.,ord1,phi1,dphi1);//g e dg             0<=x1<=1 -> -1<=2*x1-1<=1
		TPZShapeLinear::fOrthogonal(x[2],ord2,phi2,dphi2);//h e dh      -1<=x3<=1
		int index = 0;//x � ponto de integra��o dentro da pir�mide
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
    
    void TPZShapePrism::ShapeInternal(int side, TPZVec<REAL> &x, int order,
                                      TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){
        if (side < 6 || side > 20) {
            DebugStop();
        }
        
        switch (side) {
                
            case 6:
            case 7:
            case 8:
            case 9:
            case 10:
            case 11:
            case 12:
            case 13:
            case 14:
            {
                pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
            }
                break;
                
            case 15:
            case 19:
                
            {
            
                pzshape::TPZShapeTriang::ShapeInternal(x, order, phi, dphi);
            }
                break;
                
            case 16:
            case 17:
            case 18:
            {
                pzshape::TPZShapeQuad::ShapeInternal(x, order, phi, dphi);
            }
                break;

            case 20:
            {
                ShapeInternal(x, order, phi, dphi);
            }
                break;
            default:
                std::cout << "Wrong side parameter side " << side << std::endl;
                DebugStop();
                break;
        }
     
        
    }
	
	void TPZShapePrism::TransformDerivativeFromRibToPrisma(int rib,int num,TPZFMatrix<REAL> &dphi) {
		for (int j = 0;j<num;j++) {
			dphi(2,j) = gRibTrans3dPrisma1d[rib][2]*dphi(0,j);
			dphi(1,j) = gRibTrans3dPrisma1d[rib][1]*dphi(0,j);
			dphi(0,j) = gRibTrans3dPrisma1d[rib][0]*dphi(0,j);
		}
	}
	
	void TPZShapePrism::TransformDerivativeFromFaceToPrisma(int face,int num,TPZFMatrix<REAL> &dphi) {
		for (int j = 0;j<num;j++) {
			dphi(2,j) = gFaceTrans3dPrisma2d[face][0][2]*dphi(0,j)+gFaceTrans3dPrisma2d[face][1][2]*dphi(1,j);
			REAL dphi1j = dphi(1,j);
			dphi(1,j) = gFaceTrans3dPrisma2d[face][0][1]*dphi(0,j)+gFaceTrans3dPrisma2d[face][1][1]*dphi(1,j);
			dphi(0,j) = gFaceTrans3dPrisma2d[face][0][0]*dphi(0,j)+gFaceTrans3dPrisma2d[face][1][0]*dphi1j;//dphi(1,j);
		}
	}
	//transforma a derivada no ponto dentro da face
	void TPZShapePrism::TransformDerivativeFace3dPrisma(int transid, int face, int num, TPZFMatrix<REAL> &in) {
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
	
	int TPZShapePrism::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}

};
