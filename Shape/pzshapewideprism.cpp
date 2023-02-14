/**
 * @file
 * @brief Contains the implementation of the TPZShapeWidePrism methods.
 */

#include "pzshapewideprism.h"
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
	
void TPZShapeWidePrism::InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders)
{
    int nshape = (order-1)*(order)*(order-1)/2;
    if (shapeorders.Rows() != nshape) {
        DebugStop();
    }
    int count = 0;
    int ord1 = order - 1;
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
}

	
	int TPZShapeWidePrism::NConnectShapeF(int side, int order) {
		if(side<6) return 1;//0 a 4
		if(side<15) return (order-1);//6 a 14
		if(side==15 || side==19) {
			return ((order-2)*(order-1)/2);
		}
		if(side>15 && side<19) {//16,17,18
			return ((order-1)*(order-1));
		}
		if(side==20) {
			return ((order-1)*(order)*(order-1)/2);
		}
		PZError << "TPZShapeWidePrism::NConnectShapeF, bad parameter side " << side << endl;
		return 0;
	}
	
	int TPZShapeWidePrism::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}

};
