/**
 * @file
 * @brief Contains the implementation of the TPZShapePiram methods.
 */

#include "pzshapepiramHdiv.h"


using namespace std;

namespace pzshape {
    
    void TPZShapePiramHdiv::Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
        
        // not supported anymore
        DebugStop();
    }
    
    int TPZShapePiramHdiv::NConnectShapeF(int side, int order) {
        switch (side) {
            case 0:
            case 1:
            case 2:
            case 3:
            case 5:
            case 6:
            case 7:
            case 8:
            case 13:
                
                return TPZShapePiram::NConnectShapeF(side, order)*2;
                break;
            case 4:
                return 0;
                break;
            default:
                return TPZShapePiram::NConnectShapeF(side, order);
                break;
        }
    }
    
    
};
