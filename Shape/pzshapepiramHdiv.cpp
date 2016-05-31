/**
 * @file
 * @brief Contains the implementation of the TPZShapePiram methods.
 */

#include "pzshapepiramHdiv.h"


using namespace std;

namespace pzshape {
	
    void TPZShapePiramHdiv::Shape(TPZVec<REAL> &pt, TPZVec<long> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {

        REAL temp = 0.;
        if (IsZero(pt[2] - 1.)){
            temp = pt[2];
            pt[2] -= 1.e-8;
        }
        
        TPZShapePiram::Shape(pt, id, order, phi, dphi);

        // we need to push more functions, according to the order
        int allorder = order[0];
        int numshift = (allorder+1)*(allorder+1);
        int numbubbles = (allorder-1)*(allorder-1);
        
        TPZFNMatrix<20,REAL> phicp(phi);
        TPZFNMatrix<60,REAL> dphicp(dphi);
        // remove the tip function
        for (int i = 4; i > 0; i--) {
            phicp(i) = phicp(i-1);
            for (int j=0; j<3; j++) {
                dphicp(j,i) = dphicp(j,i-1);
            }
        }
        
        DebugStop();
        
        // I decided to put the four new functions (without * (1-z)) in the begining of the vector
        const int stride = 4;
        for (int i = 0 ; i < phicp.Rows() - 4; i++) {
            phi(i+stride,0) = phicp(i+1,0);
            for (int j = 0; j < 3; j++) {
                dphi(j,i+stride) = dphicp(j,i+1);
            }
        }
        
        for (int i = 0; i < 4; i++) {
            phi(i,0) = phi(i+stride,0)/(1.-pt[2]);
            for (int j = 0; j < 3; j++) {
                dphi(j,i) = dphi(j,i+stride)/(1.-pt[2]);
            }
            dphi(2,i) += phi(i,0)/(1-pt[2])/(1.-pt[2]);
        }
        if (!IsZero(temp)){
            pt[2] = temp;
        }
    }
	
    int TPZShapePiramHdiv::NConnectShapeF(int side, int order) {
        if (side == 0){ // the first connect has now 4 shapefunctions
            return 4;
        } // the other connects remain the same
        else
        {
            return TPZShapePiram::NConnectShapeF(side, order);
        }
    }

    
};
