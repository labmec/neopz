#ifndef TPZReservoirDATAH
#define TPZReservoirDATAH
/*
 *  ReservoirData.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

class ReservoirData {
    
public:
    
    /** @brief State: Stiffness or Mass Matrix Calculations */
    enum EState { ELastState = 0, ECurrentState = 1 };
    
    /**
     * @ingroup Characteristic Parameters
     * @brief Define characteristic parameters for Darcy linear flow.
     * @since December 08, 2014
     */
    
    /** @brief Characteristic length - m */
    REAL fLref;
    
    /** @brief Characteristic Permeability - m2 */
    REAL fKref;
    
    /** @brief Characteristic Pressure - Pa */
    REAL fPref;
    
    /** @brief Characteristic Density - kg/m3 */
    REAL fRhoref;
    
    /** @brief Characteristic viscosity - Pa s */
    REAL fEtaref;
    
    /** @brief Layer average thickness - m */
    REAL fh;
    
    /** @brief Layer average radius - m */
    REAL fr;
    
    /** @brief well radius - m */
    REAL frw;
    
    /** @brief Layer Top depth  - m */
    REAL fDepth;
    
    /** @brief Porosity at P of reference - */
    REAL fPhiref;
    
    /** @brief Rock Compressibility 1/pa - */
    REAL fcrock;
    
    /** @brief Fluid Compressibility 1/pa - */
    REAL fcfluid;
    
    /** @brief Is GID geometry - */
    bool fIsGIDGeometry;
    
    /** @brief absolute permeability */
    TPZFMatrix<REAL> fKab;
    
    /** @brief absolute permeability inverse */
    TPZFMatrix<REAL> fKabinv;
    
    ReservoirData();
    
    ~ReservoirData();
    
    /**
     * @brief \f$ Rock porosity. \f$ Phi = Phi( P ) \f$
     * @param P fluid pressure
     */
    void Porosity(REAL P, REAL &poros, REAL &dPorosDp);
    
    /**
     * @brief \f$ Oil density RhoOil = RhoOil( P ) \f$
     * @param P fluid pressure
     */
    void Density(REAL P, REAL &Rho, REAL &dRhoDpo);
    
    /**
     * @brief Oil viscosity. \f$ OilViscosity = ViscOil( P ) \f$
     * @param P fluid pressure
     */
    void Viscosity(REAL P, REAL &Viscosity, REAL &dViscosityDpo);
    
    /** @brief Set the characteristic length - m */
    void SetLref(REAL Lref) {fLref = Lref; }
    
    /** @brief Characteristic length - m */
    REAL Lref() {return fLref; }
    
    /** @brief Set the geometry source */
    void SetIsGIDGeometry(bool isGIDgeom) {fIsGIDGeometry = isGIDgeom; }
    
    /** @brief Get the geometry source */
    bool GetIsGIDGeometry() {return fIsGIDGeometry; }
    
    /** @brief Set the average thickness - m */
    void SetLayerTop(REAL D) {fDepth = D; }
    
    /** @brief Get the average thickness - m */
    REAL LayerTop() {return fDepth; }
    
    /** @brief Set the average thickness - m */
    void SetLayerh(REAL h) {fh = h; }
    
    /** @brief Get the average thickness - m */
    REAL Layerh() {return fh; }
    
    /** @brief Set the average radius - m */
    void SetLayerr(REAL r) {fr = r; }
    
    /** @brief Get the average radius - m */
    REAL Layerr() {return fr; }
    
    /** @brief Set the well radius - m */
    void SetLayerrw(REAL rw) {frw = rw; }
    
    /** @brief Get the well radius - m */
    REAL Layerrw() {return frw; }
    
    /** @brief Set the characteristic Permeability - m2 */
    void SetKref(REAL Kref) {fKref = Kref;}
    
    /** @brief Characteristic Permeability - m2 */
    REAL Kref() {return fKref;}
    
    /** @brief Set the characteristic Pressure - Pa */
    void SetPref(REAL Pref) {fPref = Pref;}
    
    /** @brief Characteristic Pressure - Pa */
    REAL Pref() {return fPref;}
    
    /** @brief Set the characteristic Density - kg/m3 */
    void Rhoref(REAL Rhoref) {fRhoref = Rhoref;}
    
    /** @brief Characteristic Density - kg/m3 */
    REAL Rhoref() {return fRhoref;}
    
    /** @brief Set the characteristic viscosity - Pa s */
    void SetEtaref(REAL Etaref) {fEtaref = Etaref;}
    
    /** @brief Characteristic viscosity - Pa s */
    REAL Etaref() {return fEtaref;}
    
    /** @brief Porosity at P of reference - */
    void SetPhiRef(REAL Phiref) {fPhiref = Phiref;}
    
    /** @brief Porosity at P of reference - */
    REAL PhiRef() {return fPhiref;}
    
    /** @brief Set rock compressibility 1/pa- */
    void SetcRock(REAL cr) {fcrock = cr;}
    
    /** @brief Get rock compressibility 1/pa - */
    REAL CRock() {return fcrock;}
    
    /** @brief Set fluid compressibility 1/pa- */
    void SetcFluid(REAL cf) {fcfluid = cf;}
    
    /** @brief Get fluid compressibility 1/pa - */
    REAL CFluid() {return fcfluid;}
    
    /** @brief Material indexes */
    TPZVec<int> fmaterialIds;
    
    /** @brief Set the absolute Permeability - m2 */
    void SetKabsolute(TPZFMatrix<REAL> &Kab)
    {
        fKab = Kab;
        fKabinv = Kab;
        fKabinv.Zero();
        STATE detKab;
        detKab = fKab(0,0)*fKab(1,1)-fKab(1,0)*fKab(0,1);
        
        if (fabs(detKab) <= 1.0*10-14) {
            std::cout << " Kab Matrix doesn't have Inverse, det =  " << detKab << std::endl;
            DebugStop();
        }
        
        fKabinv(0,0) = +1.0*fKab(1,1)/detKab;
        fKabinv(0,1) = -1.0*fKab(0,1)/detKab;
        fKabinv(1,0) = -1.0*fKab(1,0)/detKab;
        fKabinv(1,1) = +1.0*fKab(0,0)/detKab;
        
    }
    
    /** @brief Get the absolute Permeability - m2 */
    TPZFMatrix<REAL> Kabsolute() {return fKab;}
    
    /** @brief Get the absolute Permeability inverse - 1/m2 */
    TPZFMatrix<REAL> KabsoluteInv() {return fKabinv;}
    
    /** @brief Get the material indexes */
    void SetMatIDs(TPZVec<int> &matids) {fmaterialIds=matids;}
    
    /** @brief Get the material indexes */
    TPZVec<int> GetMatIDs() {return fmaterialIds;}
    
};


#endif