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
#include <array>
#include <algorithm>

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
        TPZVec<REAL> cart(3);
        FromPrincipalToHWCart(PrincipalCoords, cart);
        HWCyl[0] = cart[0];
        HWCyl[1] = sqrt(cart[1] * cart[1] + cart[2] * cart[2]);
        HWCyl[2] = atan2(cart[2], cart[1]);
    }

    /**
     * Transforms from eigenvalues (sigma_1, sigma_2, sigma_3) to HW cartesian coordinates (sigma_1^*, sigma_2^*, sigma_3^*)
     * @param PrincipalCoords
     * @param HWCart
     */
    static void FromPrincipalToHWCart(const TPZVec<REAL> &PrincipalCoords, TPZVec<REAL> &HWCart) {
        TPZFNMatrix<9, STATE> Rot(3, 3, 0.), temp(3, 1, 0.), cart(3, 1, 0.);
        HWCart.Resize(3, 0.);
        std::array<REAL, 3> temp_array = { PrincipalCoords[0], PrincipalCoords[1], PrincipalCoords[2] };
        std::sort(temp_array.begin(), temp_array.end());
        temp(0, 0) = temp_array[2];
        temp(1, 0) = temp_array[1];
        temp(2, 0) = temp_array[0];
        GetRotMatrix(Rot);
        Rot.Multiply(temp, cart);
        HWCart[0] = cart(0, 0);
        HWCart[1] = cart(1, 0);
        HWCart[2] = cart(2, 0);
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

private:
    TPZHWTools();
    TPZHWTools(const TPZHWTools& orig);
    virtual ~TPZHWTools();

    /// Computes the rotation matrix

    static void GetRotMatrix(TPZFMatrix<STATE> &Rot) {
        const STATE SQRT1_3 = 1. / sqrt(3.);
        const STATE SQRT1_6 = 1. / sqrt(6.);
        Rot.Resize(3, 3);
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
};

#endif /* TPZHAIGHWESTERGAARDTOOLS_H */

