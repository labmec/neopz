/**
 * @file
 * @brief Contains the implementation of the TPZShapePiram methods.
 */

#include "pzshapepiramHdiv.h"


using namespace std;

namespace pzshape {
    
    void TPZShapePiramHdiv::Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
        
        REAL temp = 0.;
        if (IsZero(pt[2] - 1.)){
            temp = pt[2];
            pt[2] -= 1.e-8;
        }
        
        TPZFNMatrix<20,REAL> phicp(phi);
        TPZFNMatrix<60,REAL> dphicp(dphi);
        
        TPZShapePiram::Shape(pt, id, order, phicp, dphicp);
        
        REAL zst = 1./(1.-pt[2]);
        int count = 0;
        int countorig = 0;
        for (int is=0; is<NSides; is++) {
            if (is == 4) {
                countorig++;
                continue;
            }
            int norig = 1;
            int nshape = 2;
            if (is>4)
            {
                norig = TPZShapePiram::NConnectShapeF(is, order[is-5]);
                nshape = NConnectShapeF(is, order[is-5]);
            }
            for (int shape = 0; shape<norig; shape++)
            {
                phi(shape+count,0) = phicp(shape+countorig,0);
                for (int d=0; d<3; d++) {
                    dphi(d,shape+count) = dphicp(d,shape+countorig);
                }
            }
            if(norig*2 == nshape)
            {
                for (int shape = 0; shape<norig; shape++)
                {
                    phi(shape+norig+count,0) = phicp(shape+countorig,0)*zst;
                    for (int d=0; d<3; d++) {
                        dphi(d,shape+norig+count) = dphicp(d,shape+countorig)*zst;
                    }
                    dphi(2,shape+norig+count) += phicp(shape+countorig,0)*zst*zst;
                }
            }
            countorig+= norig;
            count += nshape;
        }
        if (IsZero(temp - 1.))
        {
            pt[2] = temp;
        }
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
    
    int TPZShapePiramHdiv::ClassId() const{
        return Hash("TPZShapePiramHdiv") ^ TPZShapePiram::ClassId() << 1;
    }
	
    
};
