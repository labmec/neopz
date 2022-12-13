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

void TPZShapePiram::InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders)
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
	
};
