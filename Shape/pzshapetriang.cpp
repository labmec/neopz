/**
 * @file
 * @brief Contains the implementation of the TPZShapeTriang methods.
 */

#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

using namespace std;

namespace pzshape {
	
	/**Transformation of the point within a triangular face */
	REAL TPZShapeTriang::gTrans2dT[6][2][2] = {//s* , t*
		{ { 1., 0.},{ 0., 1.} },
		{ { 0., 1.},{ 1., 0.} },
		{ { 0., 1.},{-1.,-1.} },//s* = t   t* = -s-t-1 ,  etc
		{ {-1.,-1.},{ 0., 1.} },
		{ {-1.,-1.},{ 1., 0.} },
		{ { 1., 0.},{-1.,-1.} }
	};
	
	REAL TPZShapeTriang::gVet2dT[6][2] = {  {0.,0.},{0.,0.},{0.,1.},{1.,0.},{1.,0.},{0.,1.} };
	
	REAL TPZShapeTriang::gRibTrans2dT1d[3][2] = { {2.,1.},{-1.,1.},{-1.,-2.} };//Cedric : 06/03/99
	
	REAL TPZShapeTriang::gVet1dT[3] = {-1.,0.,1.};
	
	
	void TPZShapeTriang::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
		
		phi(0,0) =  1.-pt[0]-pt[1];
		phi(1,0) =  pt[0];
		phi(2,0) =  pt[1];
		dphi(0,0) = -1.;
		dphi(1,0) = -1.;
		dphi(0,1) =  1.;
		dphi(1,1) =  0.;
		dphi(0,2) =  0.;
		dphi(1,2) =  1.;
	}
	
	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input/output) value of the (4) shape functions
	 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
	 */
	void TPZShapeTriang::ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		int is;
		for(is=3; is<6; is++)
		{
			int is1 = is%3;
			int is2 = (is+1)%3;
			phi(is,0) = phi(is1,0)*phi(is2,0);
			dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
			dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
		}
		int is1 = 0;
		int is2 = 1;
		int is3 = 2;
		phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
		dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
		dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);

		// Make the generating shape functions linear and unitary
		REAL mult[] = {1.,1.,1.,4.,4.,4.,27.};
		for(is=3;is<NSides; is++)
		{
			phi(is,0) *= mult[is];
			dphi(0,is) *= mult[is];
			dphi(1,is) *= mult[is];
		}

		
	}
	
	void TPZShapeTriang::Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,
							   TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		ShapeCorner(pt,phi,dphi);
		if (order[0] < 2 && order[1] < 2 && order[2] < 2 && order[3] < 3) return;
		int is,d;
		
//        TPZFNMatrix<10,REAL> phiblend(NSides,1),dphiblend(Dimension,NSides);
        REAL store1[20],store2[40];
        TPZFMatrix<REAL> phiblend(NSides,1,store1,20),dphiblend(Dimension,NSides,store2,40);
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
		int shape = 3;
		for (int rib = 0; rib < 3; rib++) {
			
			ProjectPoint2dTriangToRib(rib,pt,out);
			TPZManVector<REAL,1> outvec(1,out);
			TPZVec<int64_t> ids(2);
			ids[0] = id[rib%3];
			ids[1] = id[(rib+1)%3];
			REAL store1[20],store2[40];
			int ord2 = order[rib]-1;//numero de shapes por lado rib
			TPZFMatrix<REAL> phin(ord2,1,store1,20),dphin(2,ord2,store2,40);
			TPZShapeLinear *shplin=0;
			shplin->ShapeInternal(outvec,order[rib],phin,dphin,shplin->GetTransformId1d(ids));
			TransformDerivativeFromRibToTriang(rib,ord2,dphin);
			for (int i = 0; i < ord2; i++) {
				phi(shape,0) = phiblend(rib+3,0)*phin(i,0);
				for(int xj=0;xj<2;xj++) {
					dphi(xj,shape) = dphiblend(xj,rib+3)* phin( i, 0)+
					phiblend(rib+3, 0 )* dphin(xj,i);
				}
				shape++;
			}
		}
		if (order[3] < 3) return;//ordem na face
//        REAL store1[20],store2[40];
		int ord =  order[3]-2;//num de shapes da face
		int nsh = (ord*(ord+1))/2;
        TPZFMatrix<REAL> phin(nsh,1,store1,20),dphin(2,nsh,store2,40);
		ShapeInternal(pt,order[3]-2,phin,dphin,GetTransformId2dT(id));
		for(int i=0;i<nsh;i++)	{//number of internal shape equal maximal order
			phi(shape,0) = phiblend(6,0)*phin(i,0);
			for(int d=0;d<2;d++) {
				dphi(d,shape) = dphiblend(d,6)* phin(i,0) +
				phiblend(6,0)*dphin(d,i);
			}
			shape++;
		}
	}
	
	void TPZShapeTriang::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,
								   TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		if(side<0 || side>6) PZError << "TPZShapeTriang::SideShape. Bad paramenter side.\n";
		else if(side==6) Shape(pt,id,order,phi,dphi);
		else if(side<3) {
			TPZShapePoint::Shape(pt,id,order,phi,dphi);
		} else {
			TPZShapeLinear::Shape(pt,id,order,phi,dphi);
		}
		
		
	}
    
    void TPZShapeTriang::ShapeOrder(TPZVec<int64_t> &id, TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        int64_t nsides = TPZTriangle::NSides;
        // o que eh o vetor order?
        // Eu suponho que em cada posicao tem a ordem de cada lado.
        // Na shape ja esta associado a lados com dimensao maior que 1, order[0](lado 3) ...

        int nshape = NCornerNodes;
        
        for (int is = NCornerNodes; is < nsides; is++)
        {
            nshape += NConnectShapeF(is,order[is-NCornerNodes]);
        }
        // shapeorders ja vem com tamanho definido, entao a conta acima so serve para ver se as ordens estao corretas
        if (nshape!=shapeorders.Rows())
        {
            // Algo errado
            DebugStop();
        }
        
        int linha = 0;
        for (int side = 0; side < nsides; side++)
        {

            //nshape = NConnectShapeF(side,order[side]);
            //int sideorder = order[side];
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
    
    
    void TPZShapeTriang::SideShapeOrder(int side,  TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders)
    {
        
        if (side<=2)
        {
            if(shapeorders.Rows() != 1)
            {
                DebugStop();
            }
            shapeorders(0,0) = 1;
            shapeorders(0,1) = 0;
            shapeorders(0,2) = 0;
        }
        else if (side == 6)
        {
            int nshape = (order-2)*(order-1)/2;
            if(shapeorders.Rows() != nshape) DebugStop();
            int nsh=1;
            int ish = 0;
            for (int i=3; i<=order; i++) {
                for (int j=0; j<nsh; j++) {
                    shapeorders(ish,0) = i;
                    shapeorders(ish,1) = i;
                    ish++;
                }
                nsh++;
            }
        }
        else 
        {
            int nshape = order-1;
            if (shapeorders.Rows() != nshape) {
                DebugStop();
            }
            for (int ioy = 0; ioy < order-1; ioy++)
            {
                shapeorders(ioy,0) = ioy+2;
            }
        }
        
    }
    
    
	void TPZShapeTriang::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
									   TPZFMatrix<REAL> &dphi,int triangle_transformation_index) {
		
		if(order < 0) return;
		int ord1 = order;
		int numshape = (ord1*(ord1+1))/2;
		TPZManVector<REAL,2> out(2);
		TransformPoint2dT(triangle_transformation_index,x,out);
        
        out[0] = 2.*out[0]-1.;
        out[1] = 2.*out[1]-1.;
		
		if (phi.Rows() < numshape || dphi.Cols() < numshape) {
			PZError << "\nTPZCompEl::Shape2dTriangleInternal phi or dphi resized\n";
			phi.Resize(numshape,1);
			dphi.Resize(dphi.Rows(),numshape);
		}
		REAL store1[20],store2[20],store3[20],store4[20];
		TPZFMatrix<REAL> phi0(numshape,1,store1,20),phi1(numshape,1,store2,20),dphi0(1,numshape,store3,20),dphi1(1,numshape,store4,20);
		
		TPZShapeLinear::fOrthogonal(out[0],numshape,phi0,dphi0);
		TPZShapeLinear::fOrthogonal(out[1],numshape,phi1,dphi1);
		int index = 0;
		int i;
		for (int iplusj=0;iplusj<ord1;iplusj++) {
			for (int j=0;j<=iplusj;j++) {
				i = iplusj-j;
				phi(index,0) = phi0(i,0)*phi1(j,0);
				dphi(0,index) = 2.*dphi0(0,i)*phi1(j,0);
				dphi(1,index) = 2.*phi0(i,0)*dphi1(0,j);
				index++;
			}
		}
		TransformDerivative2dT(triangle_transformation_index,numshape,dphi);
	}
	
	void TPZShapeTriang::ProjectPoint2dTriangToRib(int rib, TPZVec<REAL> &in, REAL &out) {
		
		out = gRibTrans2dT1d[rib][0]*in[0]+gRibTrans2dT1d[rib][1]*in[1]+gVet1dT[rib];
	}
	
	void TPZShapeTriang::TransformDerivativeFromRibToTriang(int rib,int num,TPZFMatrix<REAL> &dphi) {
		
		for (int j = 0;j<num;j++) {
			
			dphi(1,j) = gRibTrans2dT1d[rib][1]*dphi(0,j);
			dphi(0,j) = gRibTrans2dT1d[rib][0]*dphi(0,j);
		}
	}
	
	int TPZShapeTriang::GetTransformId2dT(TPZVec<int64_t> &id) {
		
		int id0,id1,minid;
		id0 = (id[0] < id[1]) ? 0 : 1;
		minid = (id[2] < id[id0]) ? 2 : id0;
		id0 = (minid+1)%3;
		id1 = (minid+2)%3;
		
		if (id[id0] < id[id1]) {//antihorario
			
			if (minid == 0) return 0;
			if (minid == 1) return 2;
			if (minid == 2) return 4;
			
		} else {//horario
			
			if (minid == 0) return 1;
			if (minid == 1) return 3;
			if (minid == 2) return 5;
		}
		return 0;
	}
	
	//transf. o ponto dentro da face triangular
	void TPZShapeTriang::TransformPoint2dT(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out) {
		
		out[0] = gTrans2dT[transid][0][0]*in[0]+gTrans2dT[transid][0][1]*in[1]+gVet2dT[transid][0];
		out[1] = gTrans2dT[transid][1][0]*in[0]+gTrans2dT[transid][1][1]*in[1]+gVet2dT[transid][1];
	}
	
	void TPZShapeTriang::TransformDerivative2dT(int transid, int num, TPZFMatrix<REAL> &in) {
		
		int i;
		for(i=0;i<num;i++) { //ds/dcsi
			REAL aux[2];
			aux[0] = in(0,i);
			aux[1] = in(1,i);
			in(0,i) = gTrans2dT[transid][0][0]*aux[0]+gTrans2dT[transid][1][0]*aux[1];
			in(1,i) = gTrans2dT[transid][0][1]*aux[0]+gTrans2dT[transid][1][1]*aux[1];
		}
	}
	
	int TPZShapeTriang::NConnectShapeF(int side, int order) {
		switch(side) {
			case 0:
			case 1:
			case 2:
				return 1;
			case 3:
			case 4:
			case 5:
				return order-1;
			case 6:
				return (order-2) < 0 ? 0 : ((order-2)*(order-1))/2;
			default:
				PZError << "TPZShapeTriang::NConnectShapeF, bad parameter iconnect " << side << endl;
				return 0;
		}
	}
	
	int TPZShapeTriang::NShapeF(TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++)
        {
            res += NConnectShapeF(in,order[in-NCornerNodes]);
        }
		return res;
	}

};
