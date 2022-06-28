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
	

	/*
	void TPZShapeCube::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
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
    */
    void TPZShapeCube::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        //DebugStop();
    
        int64_t nsides = TPZShapeCube::NSides;
        int nshape = NCornerNodes;
        
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
    
    
    void TPZShapeCube::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        //DebugStop();
        if (side<=7)
        {
            if (shapeorders.Rows() != 1)
            {
                DebugStop();
            }
            shapeorders(0,0) = 1;
            shapeorders(0,1) = 0;
            shapeorders(0,2) = 0;
        }
        else if (side == 26)
        {
            int nsei = (order-1)*(order-1)*(order-1);
            if (shapeorders.Rows() != nsei) {
                DebugStop();
            }
            int count = 0;
            for (int i=2; i<order+1; i++) {
                for (int j=2; j<order+1; j++) {
                    for (int k=2; k<order+1; k++) {
                        int a = i;
                        int b = j;
                        int c = k;
                        shapeorders(count,0) = a;
                        shapeorders(count,1) = b;
                        shapeorders(count,2) = c;
                        count++;

                    }
                    
                }
            }
        }
        else if (side>7 && side<20)
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
        else
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
	
	int TPZShapeCube::NShapeF(const TPZVec<int> &order) {
		int in,res = NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	

};
