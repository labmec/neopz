/**
 * @file
 * @brief Contains the implementation of the TPZShapeQuad methods.
 */

#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.shape.TPZShapeQuad");
#endif

using namespace std;

namespace pzshape {
	
	
    void TPZShapeQuad::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        int64_t nsides = TPZQuadrilateral::NSides;
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
    
void TPZShapeQuad::InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders)
{
    if (shapeorders.Rows() != (order-1)*(order-1)) {
        DebugStop();
    }
    int transid = GetTransformId(id);
    int count = 0;
    for (int i=2; i<order+1; i++) {
        for (int j=2; j<order+1; j++) {
            int a = i;
            int b = j;
            switch (transid)
            {
                case 0:
                case 3:
                case 4:
                case 7:
                    break;
                case 1:
                case 2:
                case 5:
                case 6:
                {
                    int c = a;
                    a = b;
                    b = c;
                }
                    break;
                    
                default:
                    DebugStop();
                    break;
            }
            shapeorders(count,0) = a;
            shapeorders(count,1) = b;
            count++;
        }
    }
}

    void TPZShapeQuad::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        
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
        else if (side == 8)
        {
            InternalShapeOrder(id, order, shapeorders);
        }
        else
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


    }
    
	
int TPZShapeQuad::NConnectShapeF(int side, int order) {
    if(order <1) order++;
#if PZDEBUG
    if(order < 1){
      PZError << "TPZShapeCube::NConnectShapeF, bad parameter order " << order << endl;
      DebugStop();
    }
#endif
		if(side<4) return 1;//0 a 4
		//   int s = side-4;//s = 0 a 14 ou side = 6 a 20
#ifdef PZDEBUG
        {
//            if(order <1) DebugStop();
        }
#endif
		if(side<8) return (order-1);//6 a 14
		if(side==8) {
			return ((order-1)*(order-1));
		}
		PZError << "TPZShapeQuad::NConnectShapeF, bad parameter side " << side << endl;
    DebugStop();
		return 0;
	}
	
	int TPZShapeQuad::NShapeF(const TPZVec<int> &order) {
			
				
		int in,res=NCornerNodes;
			for(in=NCornerNodes;in<NSides;in++) {
									
					res += NConnectShapeF(in,order[in-NCornerNodes]);
			}
		return res;
	}
	
	
	

};



