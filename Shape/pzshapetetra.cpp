/**
 * @file
 * @brief Contains the implementation of the TPZShapeTetra methods.
 */

#include "pzshapetetra.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

using namespace std;

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
	
    void TPZShapeTetra::ShapeCorner(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
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
    
    
    //troco para ShapeCorner
	void TPZShapeTetra::CornerShape(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
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
	
	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input/output) value of the (4) shape functions
	 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
	 */
	void TPZShapeTetra::ShapeGenerating(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		int is;
		// 6 ribs
		for(is=4; is<NSides; is++)
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
					//int face = is-10;
					int is1 = SideNodeLocId(is,0); //ShapeFaceId[face][0]; 
					int is2 = SideNodeLocId(is,1); //ShapeFaceId[face][1]; 
					int is3 = SideNodeLocId(is,2); //ShapeFaceId[face][2]; 
					phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
					dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
					dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
					dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
				}
					break;
				case 4:
				{
					phi(is,0) = phi(0,0)*phi(1,0)*phi(2,0)*phi(3,0);
					for(int xj=0;xj<3;xj++) {
						dphi(xj,is) = dphi(xj,0)* phi(1 ,0)* phi(2 ,0)* phi(3 ,0) +
						phi(0, 0)*dphi(xj,1)* phi(2 ,0)* phi(3 ,0) +
						phi(0, 0)* phi(1 ,0)*dphi(xj,2)* phi(3 ,0) +
						phi(0, 0)* phi(1 ,0)* phi(2 ,0)*dphi(xj,3);
					}
				}
					break;
					
				default:
					DebugStop();
			}
		}

		REAL mult[] = {1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,27.,27.,27.,27.,54.};
		for(is=4;is<NSides; is++)
		{
			phi(is,0) *= mult[is];
			dphi(0,is) *= mult[is];
			dphi(1,is) *= mult[is];
			dphi(2,is) *= mult[is];
		}
		 
	}
	
/**
 * @brief Computes the generating shape functions for a quadrilateral element
 * @param pt (input) point where the shape function is computed
 * @param phi (input/output) value of the (4) shape functions
 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
 */
void TPZShapeTetra::ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
    REAL mult[] = {1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,27.,27.,27.,27.,54.};

    // 6 ribs
    for(int is=4; is<NSides; is++)
    {
        if(nshape[is-4] < 1) continue;
        int nsnodes = NSideNodes(is);
        switch(nsnodes)
        {
            case 2:
            {
                int is1 = SideNodeLocId(is,0);
                int is2 = SideNodeLocId(is,1);
                phi(is,0) = mult[is]*phi(is1,0)*phi(is2,0);
                dphi(0,is) = mult[is]*(dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2));
                dphi(1,is) = mult[is]*(dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2));
                dphi(2,is) = mult[is]*(dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2));
            }
                break;
            case 3:
            {
                //int face = is-10;
                int is1 = SideNodeLocId(is,0); //ShapeFaceId[face][0];
                int is2 = SideNodeLocId(is,1); //ShapeFaceId[face][1];
                int is3 = SideNodeLocId(is,2); //ShapeFaceId[face][2];
                phi(is,0) = mult[is]*phi(is1,0)*phi(is2,0)*phi(is3,0);
                dphi(0,is) = mult[is]*(dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3));
                dphi(1,is) = mult[is]*(dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3));
                dphi(2,is) = mult[is]*(dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3));
            }
                break;
            case 4:
            {
                phi(is,0) = mult[is]*phi(0,0)*phi(1,0)*phi(2,0)*phi(3,0);
                for(int xj=0;xj<3;xj++) {
                    dphi(xj,is) = mult[is]*(dphi(xj,0)* phi(1 ,0)* phi(2 ,0)* phi(3 ,0) +
                    phi(0, 0)*dphi(xj,1)* phi(2 ,0)* phi(3 ,0) +
                    phi(0, 0)* phi(1 ,0)*dphi(xj,2)* phi(3 ,0) +
                    phi(0, 0)* phi(1 ,0)* phi(2 ,0)*dphi(xj,3));
                }
            }
                break;
                
            default:
                DebugStop();
        }
    }     
}


	
	//ifstream inn("mats.dt");
	void TPZShapeTetra::Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
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
		
		//  if(order[9]<2) return;
		int shape = 4;
		//rib shapes
		for (int rib = 0; rib < 6; rib++) {
			REAL outval;
			ProjectPoint3dTetraToRib(rib,pt,outval);
			TPZVec<int64_t> ids(2);
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
			TransformDerivativeFromRibToTetra(rib,ordin,dphin);
			for (int i = 0; i < ordin; i++) {
				phi(shape,0) = phiblend(rib+4,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj ,rib+4) * phin( i, 0) +
					phiblend(rib+4, 0 )  * dphin(xj,i);
				}
				shape++;
			}
		}
		//  if(order[10]<3) return;
		//face shapes
		for (int face = 0; face < 4; face++) {
			if (order[6+face] < 3) continue;
			TPZManVector<REAL,2> outval(2);
			ProjectPoint3dTetraToFace(face,pt,outval);
			REAL store1[20],store2[60];
			int ord1;//,ord2;
			//elt->FaceOrder(face,ord1,ord2);
			ord1 = order[6+face];
			//ord2 = ord1;
			if(ord1<3) continue;
			int ordin =  (ord1-2)*(ord1-1)/2;
			TPZFMatrix<REAL> phin(ordin,1,store1,20),dphin(3,ordin,store2,60);//ponto na face
			phin.Zero();
			dphin.Zero();
			TPZManVector<int64_t> ids(3);
			//	int id0,id1,id2;
			int i;
			for(i=0;i<3;i++) ids[i] = id[FaceNodes[face][i]];
			//    id0 = ShapeFaceId[face][0];//indice das shapes da face que compoem a shape atual
			//    id1 = ShapeFaceId[face][1];//equivale a FaceIdsCube(face,ids,id,id0,id1);
			//    id2 = ShapeFaceId[face][2];
			int transid = TPZShapeTriang::GetTransformId2dT(ids);
			outval[0] = (outval[0]+1.)/2.;//devido a corre��o na fun��o
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
				phi(shape,0) = phiblend(face+10,0)*phin(i,0);
				for(int xj=0;xj<3;xj++) {
					dphi(xj,shape) = dphiblend(xj,face+10) * phin(i ,0) +
					phiblend(face+10, 0) * dphin(xj,i);
				}
				shape++;
			}
		}
		if(order[10]<4) return;
		//volume shapes
		int totsum = 0,sum;
		int i;
		for(i=0;i<order[10]-3;i++) {
			sum = (i+1)*(i+2) / 2;
			totsum += sum;
		}
		int ord = totsum;
		TPZFNMatrix<80> phin(ord,1),dphin(3,ord);
		phin.Zero();
		dphin.Zero();
		ShapeInternal(pt,order[10],phin,dphin);
		for(i=0;i<ord;i++)	{
			phi(shape,0) = phiblend(NSides-1,0)*phin(i,0);
			for(int xj=0;xj<3;xj++) {
				dphi(xj,shape) = dphiblend(xj,NSides-1) * phin(i ,0) +
				phiblend(NSides-1, 0) * dphin(xj,i);
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
	
	void TPZShapeTetra::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
									  TPZFMatrix<REAL> &dphi) {
		if(order < 4) return;
		int ord = order-3;
		
		TPZFNMatrix<100, REAL> phi0(ord,1),phi1(ord,1),phi2(ord,1),dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
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
    
    void TPZShapeTetra::ShapeInternal(int side, TPZVec<REAL> &x, int order,
                                      TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){
        if (side < 4 || side > 14) {
            DebugStop();
        }
        
        switch (side) {
                
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:
            {
                pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
            }
                break;
                
            case 10:
            case 11:
            case 12:
            case 13:
                
            {
                pzshape::TPZShapeTriang::ShapeInternal(x, order, phi, dphi);
            }
                break;
                
            case 14:
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
	
	void TPZShapeTetra::TransformDerivativeFromRibToTetra(int rib,int num,TPZFMatrix<REAL> &dphi) {
		for (int j = 0;j<num;j++) {
			dphi(2,j) = gRibTrans3dTetra1d[rib][2]*dphi(0,j);
			dphi(1,j) = gRibTrans3dTetra1d[rib][1]*dphi(0,j);
			dphi(0,j) = gRibTrans3dTetra1d[rib][0]*dphi(0,j);
		}
	}
	
	void TPZShapeTetra::TransformDerivativeFromFaceToTetra(int face,int num,TPZFMatrix<REAL> &dphi) {
		
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
	
	int TPZShapeTetra::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	
	void TPZShapeTetra::SideShape(int side, TPZVec<REAL> &point, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
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
    
    void TPZShapeTetra::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        //DebugStop();
        
        int64_t nsides = TPZShapeTetra::NSides;
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
    
    
    void TPZShapeTetra::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        //DebugStop();
        if (side<=3)
        {
            if (shapeorders.Rows() != 1)
            {
                DebugStop();
            }
            shapeorders(0,0) = 1;
            shapeorders(0,1) = 0;
            shapeorders(0,2) = 0;
        }
        else if (side>=4 && side<=9)
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
        else if (side >= 10 && side <= 13)
        {
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
        {   // interno
            int totsum = 0,sum;
            int i;
            for(i=0;i<order-3;i++) {
                sum = (i+1)*(i+2) / 2;
                totsum += sum;
            }
            int nshape = totsum;
            
            if (shapeorders.Rows() != nshape) {
                DebugStop();
            }
            int count = 0;
            int ord = order-3;
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        int a = i;
                        int b = j;
                        int c = k;
                        int soma = a+b+c;
                        if( soma < ord ) { // Duvida
                            shapeorders(count,0) = 4 + soma;
                            shapeorders(count,1) = 4 + soma;
                            shapeorders(count,2) = 4 + soma;
                            count++;
                        }
                    }
                }
            }
//            int orderplus1 = order+1;
//            for (int i=4; i<orderplus1; i++)
//            {
//                for (int j=4; j<orderplus1; j++)
//                {
//                    for (int k=4; k<orderplus1; k++)
//                    {
//                        int a = i;
//                        int b = j;
//                        int c = k;
//                        int soma = a+b+c;
//                        if ((soma)<orderplus1)// Duvida
//                        {
//                            shapeorders(count,0) = soma;
//                            shapeorders(count,1) = soma;
//                            shapeorders(count,2) = soma;
//                            count++;
//                        }
//                    }
//                    
//                }
//            }

        }
    }
	
};
