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
	
    void TPZShapeTriang::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
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

void TPZShapeTriang::InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders)
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

    
    void TPZShapeTriang::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
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
            InternalShapeOrder(id, order, shapeorders);
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
	
	int TPZShapeTriang::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++)
        {
            res += NConnectShapeF(in,order[in-NCornerNodes]);
        }
		return res;
	}

};
