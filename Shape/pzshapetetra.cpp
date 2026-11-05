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
	
	
	
	int TPZShapeTetra::NConnectShapeF(int side, int order) {
#if PZDEBUG
    if(order < 1){
      PZError << "TPZShapeCube::NConnectShapeF, bad parameter order " << order << endl;
      DebugStop();
    }
#endif
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
		PZError << "TPZShapeTetra::NConnectShapeF, bad parameter side " << side << endl;
    DebugStop();
		return 0;
	}
	
	int TPZShapeTetra::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
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
    
void TPZShapeTetra::InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders)
{
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

        }
    }
	
};
