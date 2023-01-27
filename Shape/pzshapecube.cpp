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
	
void TPZShapeCube::InternalShapeOrder(const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders) {
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

/*
    void TPZShapeCube::ShapeInternal(int side, TPZVec<REAL> &x, int order,
                                       TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
        if (side < 8 || side > 26) {
            DebugStop();
        }
        
        switch (side) {
                
            case 8:
            case 9:
            case 10:
            case 11:
            case 12:
            case 13:
            case 14:
            case 15:
            case 16:
            case 17:
            case 18:
            case 19:
            {
                pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
            }
                break;
            case 20:
            case 21:
            case 22:
            case 23:
            case 24:
            case 25:
            {
                pzshape::TPZShapeQuad::ShapeInternal(x, order, phi, dphi);
            }
                break;
            case 26:
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
*/
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
