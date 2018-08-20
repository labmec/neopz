/* 
 * File:   TPZHaighWestergaardTools.h
 * Author: quinelato
 *
 * Created on August 15, 2017, 4:28 PM
 */

#ifndef TPZHWTOOLS_H
#define TPZHWTOOLS_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <math.h>


/**
 * Performs various transformations using the Haigh-Westergaard representation.
 */
class TPZHWTools {
public:

    /**
     * Transforms from Haigh Westergaard cylindrical coordinates (xi, rho, beta) to Haigh Westergaard cartesian coordinates (sigma_1^*, sigma_2^*, sigma_3^*)
     * @param HWCylCoords
     * @param HWCart
     */
    static void FromHWCylToHWCart(const TPZVec<REAL> &HWCylCoords, TPZVec<REAL> &HWCart) {
        HWCart[0] = HWCylCoords[0];
        HWCart[1] = HWCylCoords[1] * cos(HWCylCoords[2]);
        HWCart[2] = HWCylCoords[1] * sin(HWCylCoords[2]);
    }

    /**
     * Transforms from eigenvalues (sigma_1, sigma_2, sigma_3) to HW Cylindrical coordinates (xi, rho, beta)
     * @param PrincipalCoords
     * @param HWCyl
     */
    static void FromPrincipalToHWCyl(const TPZVec<REAL> &PrincipalCoords, TPZVec<REAL> &HWCyl) {
        TPZManVector<REAL,3> cart(3);
        FromPrincipalToHWCart(PrincipalCoords, cart);
        HWCyl[0] = cart[0]; // xi
        HWCyl[1] = sqrt(cart[1] * cart[1] + cart[2] * cart[2]); // rho
        HWCyl[2] = atan2(cart[2], cart[1]); // beta
    }

    /**
     * Transforms from eigenvalues (sigma_1, sigma_2, sigma_3) to HW cartesian coordinates (sigma_1^*, sigma_2^*, sigma_3^*)
     * @param PrincipalCoords
     * @param HWCart
     */
    static void FromPrincipalToHWCart(const TPZVec<REAL> &PrincipalCoords, TPZVec<REAL> &HWCart) {
        TPZFNMatrix<9, REAL> Rot(3, 3, 0.);
//        HWCart.Resize(3, 0.); @omar::time_profiling
        GetRotMatrix(Rot);
//        Rot.Multiply(temp, cart); @omar::time_profiling
        HWCart[0] = Rot(0,0) * PrincipalCoords[0] + Rot(0,1) * PrincipalCoords[1] + Rot(0,2) * PrincipalCoords[2];
        HWCart[1] = Rot(1,0) * PrincipalCoords[0] + Rot(1,1) * PrincipalCoords[1] + Rot(1,2) * PrincipalCoords[2];
        HWCart[2] = Rot(2,0) * PrincipalCoords[0] + Rot(2,1) * PrincipalCoords[1] + Rot(2,2) * PrincipalCoords[2];
    }

    /**
     * Transforms from HW Cylindrical coordinates (xi, rho, beta) to eigenvalues (sigma_1, sigma_2, sigma_3)
     * @param HWCylCoords
     * @param PrincipalCoords
     */
    static void FromHWCylToPrincipal(const TPZVec<REAL> &HWCylCoords, TPZVec<REAL> &PrincipalCoords) {
        const REAL tempWithXi = (1. / sqrt(3.)) * HWCylCoords[0];
        const REAL tempWithRho = sqrt(2. / 3.) * HWCylCoords[1];
        PrincipalCoords[0] = tempWithXi + tempWithRho * cos(HWCylCoords[2]);
        PrincipalCoords[1] = tempWithXi + tempWithRho * cos(HWCylCoords[2]-(2. * M_PI / 3.));
        PrincipalCoords[2] = tempWithXi + tempWithRho * cos(HWCylCoords[2]+(2. * M_PI / 3.));
    }
    
    /**
     * Computes a 2x2 matrix inverse analitically being used during 2 parameters optimization
     */
    static void A2x2Inverse(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &Ainv){
        
#ifdef PZDEBUG
        if ((A.Rows() != 2 && A.Cols() != 2) || (Ainv.Rows() != 2 && Ainv.Cols() != 2)) {
            DebugStop();
        }
#endif
        Ainv(0,0)= A(1,1)/(-(A(0,1)*A(1,0)) + A(0,0)*A(1,1));
        Ainv(0,1)= -(A(0,1)/(-(A(0,1)*A(1,0)) + A(0,0)*A(1,1)));
        Ainv(1,0)= -(A(1,0)/(-(A(0,1)*A(1,0)) + A(0,0)*A(1,1)));
        Ainv(1,1)= A(0,0)/(-(A(0,1)*A(1,0)) + A(0,0)*A(1,1));
    }
    /**
     * Computes a 3x3 matrix inverse analitically being used during 3 parameters optimization
     */
    static void A3x3Inverse(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &Ainv){
        
#ifdef PZDEBUG
        if ((A.Rows() != 3 && A.Cols() != 3) || (Ainv.Rows() != 3 && Ainv.Cols() != 3)) {
            DebugStop();
        }
#endif
        Ainv(0,0)= (A(1,2)*A(2,1) - A(1,1)*A(2,2))/(A(0,2)*A(1,1)*A(2,0) - A(0,1)*A(1,2)*A(2,0) - A(0,2)*A(1,0)*A(2,1) + A(0,0)*A(1,2)*A(2,1) + A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,1)*A(2,2));
        Ainv(0,1)= (A(0,2)*A(2,1) - A(0,1)*A(2,2))/(-(A(0,2)*A(1,1)*A(2,0)) + A(0,1)*A(1,2)*A(2,0) + A(0,2)*A(1,0)*A(2,1) - A(0,0)*A(1,2)*A(2,1) - A(0,1)*A(1,0)*A(2,2) + A(0,0)*A(1,1)*A(2,2));
        Ainv(0,2)= (A(0,2)*A(1,1) - A(0,1)*A(1,2))/(A(0,2)*A(1,1)*A(2,0) - A(0,1)*A(1,2)*A(2,0) - A(0,2)*A(1,0)*A(2,1) + A(0,0)*A(1,2)*A(2,1) + A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,1)*A(2,2));
        Ainv(1,0)= (A(1,2)*A(2,0) - A(1,0)*A(2,2))/(-(A(0,2)*A(1,1)*A(2,0)) + A(0,1)*A(1,2)*A(2,0) + A(0,2)*A(1,0)*A(2,1) - A(0,0)*A(1,2)*A(2,1) - A(0,1)*A(1,0)*A(2,2) + A(0,0)*A(1,1)*A(2,2));
        Ainv(1,1)= (A(0,2)*A(2,0) - A(0,0)*A(2,2))/(A(0,2)*A(1,1)*A(2,0) - A(0,1)*A(1,2)*A(2,0) - A(0,2)*A(1,0)*A(2,1) + A(0,0)*A(1,2)*A(2,1) + A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,1)*A(2,2));
        Ainv(1,2)= (A(0,2)*A(1,0) - A(0,0)*A(1,2))/(-(A(0,2)*A(1,1)*A(2,0)) + A(0,1)*A(1,2)*A(2,0) + A(0,2)*A(1,0)*A(2,1) - A(0,0)*A(1,2)*A(2,1) - A(0,1)*A(1,0)*A(2,2) + A(0,0)*A(1,1)*A(2,2));
        Ainv(2,0)= (A(1,1)*A(2,0) - A(1,0)*A(2,1))/(A(0,2)*A(1,1)*A(2,0) - A(0,1)*A(1,2)*A(2,0) - A(0,2)*A(1,0)*A(2,1) + A(0,0)*A(1,2)*A(2,1) + A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,1)*A(2,2));
        Ainv(2,1)= (A(0,1)*A(2,0) - A(0,0)*A(2,1))/(-(A(0,2)*A(1,1)*A(2,0)) + A(0,1)*A(1,2)*A(2,0) + A(0,2)*A(1,0)*A(2,1) - A(0,0)*A(1,2)*A(2,1) - A(0,1)*A(1,0)*A(2,2) + A(0,0)*A(1,1)*A(2,2));
        Ainv(2,2)= (A(0,1)*A(1,0) - A(0,0)*A(1,1))/(A(0,2)*A(1,1)*A(2,0) - A(0,1)*A(1,2)*A(2,0) - A(0,2)*A(1,0)*A(2,1) + A(0,0)*A(1,2)*A(2,1) + A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,1)*A(2,2));
        
    }
    
    /// Computes the rotation matrix
    static void GetRotMatrix(TPZFMatrix<REAL> &Rot) {
        const REAL SQRT1_3 = 1. / sqrt(3.);
        const REAL SQRT1_6 = 1. / sqrt(6.);
        //        Rot.Resize(3, 3); // omar::time_profiling
        Rot(0, 0) = SQRT1_3;
        Rot(0, 1) = SQRT1_3;
        Rot(0, 2) = SQRT1_3;
        Rot(1, 0) = M_SQRT2*SQRT1_3;
        Rot(1, 1) = -SQRT1_6;
        Rot(1, 2) = -SQRT1_6;
        Rot(2, 0) = 0;
        Rot(2, 1) = M_SQRT1_2;
        Rot(2, 2) = -M_SQRT1_2;
    }
    
    /// Computes the inverse of the rotation matrix
    static void GetRotInvMatrix(TPZFMatrix<REAL> &RotInv) {
        TPZFMatrix<REAL> Rot(3,3);
        GetRotMatrix(Rot);
        A3x3Inverse(Rot, RotInv);
    }

public:
    
    TPZHWTools();
    
    TPZHWTools(const TPZHWTools& orig);
    
    virtual ~TPZHWTools();
    
};

#endif /* TPZHAIGHWESTERGAARDTOOLS_H */

