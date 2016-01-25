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
    
private:
    
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
    
    /** @brief Irreducible Saturation of the wetting phase - */
    REAL fS_wett_r;
    
    /** @brief Irreducible Saturation of the no wetting phase - */
    REAL fS_nwett_r;
    
    /** @brief Is GID geometry - */
    bool fIsGIDGeometry;
    
    /** @brief Definition of the initial boundaries */
    TPZStack< TPZVec<REAL> > fInitial_BCs;
    
    /** @brief Definition of the boundaries */
    TPZStack< TPZVec<REAL> > fBCs;
    
    /** @brief Absolute permeability */
    TPZFMatrix<REAL> fKab;
    
    /** @brief Absolute permeability inverse */
    TPZFMatrix<REAL> fKabinv;
    
    /** @brief Absolute permeability vector */
    TPZManVector< TPZFMatrix<REAL> > fKabVec;
   
public:
    
    ReservoirData();
    
    ~ReservoirData();
    
    /**
     * @brief \f$ Rock porosity. \f$ Phi = Phi( P ) \f$
     * @param P fluid pressure
     */
    void Porosity(REAL P, REAL &poros, REAL &dPorosDp);
    
    /** @brief Set the characteristic length - m */
    void SetLref(REAL Lref) {fLref = Lref; }
    
    /** @brief Set the characteristic Pressure - Pa */
    void SetPref(REAL Pref) {fPref = Pref;}
    
    /** @brief Get characteristic Pressure - Pa */
    REAL Pref() {return fPref;}
    
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
    
    /** @brief Set the layer thickness - m */
    void SetLayerh(REAL h) {fh = h; }
    
    /** @brief Get the layer thickness - m */
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
    
    /** @brief Get Characteristic Permeability - m2 */
    REAL Kref() {return fKref;}
    
    /** @brief Set Porosity at P of reference - */
    void SetPhiRef(REAL Phiref) {fPhiref = Phiref;}
    
    /** @brief Get Porosity at P of reference - */
    REAL PhiRef() {return fPhiref;}
    
    /** @brief Set rock compressibility 1/pa- */
    void SetcRock(REAL cr) {fcrock = cr;}
    
    /** @brief Get rock compressibility 1/pa - */
    REAL CRock() {return fcrock;}
    
    /** @brief Set the irreducible Saturation of the wetting phase - */
    void SetS_wett_r(REAL S_w_r) {fS_wett_r = S_w_r;}
    
    /** @brief Get the irreducible Saturation of the wetting phase - */
    REAL S_wett_r() {return fS_wett_r;}
    
    /** @brief Set the irreducible Saturation of the no wetting phase - */
    void SetS_nwett_r(REAL S_nw_r) {fS_nwett_r = S_nw_r;}
    
    /** @brief Get the irreducible Saturation of the no wetting phase - */
    REAL S_nwett_r() {return fS_nwett_r;}
    
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
    
    /** @brief Compute the absolute Permeability - m2 */
    TPZFMatrix<REAL> ComputeInvKabsolute(TPZFMatrix<REAL> &Kab)
    {
        TPZFMatrix<REAL> invK(2,2);
        invK.Zero();
        STATE detKab;
        detKab = Kab(0,0)*Kab(1,1)-Kab(1,0)*Kab(0,1);
        
        if (fabs(detKab) <= 1.0*10-14) {
            std::cout << " Kab Matrix doesn't have Inverse, det =  " << detKab << std::endl;
            DebugStop();
        }
        
        invK(0,0) = +1.0*Kab(1,1)/detKab;
        invK(0,1) = -1.0*Kab(0,1)/detKab;
        invK(1,0) = -1.0*Kab(1,0)/detKab;
        invK(1,1) = +1.0*Kab(0,0)/detKab;
        
        return invK;
        
    }
    
    /** @brief set Absolute permeability vector */
    void SetKvector (TPZManVector< TPZFMatrix<REAL> > KabVec) {fKabVec = KabVec; }
    
    /** @brief set Absolute permeability vector */
    TPZManVector< TPZFMatrix<REAL> > GetKvector() {return fKabVec; }
    
    /** @brief Get the absolute Permeability - m2 */
    TPZFMatrix<REAL> Kabsolute() {return fKab;}
    
    /** @brief Get the absolute Permeability inverse - 1/m2 */
    TPZFMatrix<REAL> KabsoluteInv() {return fKabinv;}

    /** @brief Get the absolute Permeability - m2 */
    TPZFMatrix<REAL> Kabsolute(TPZManVector<REAL> & x);
    
    /** @brief Get the absolute Permeability inverse - 1/m2 */
    TPZFMatrix<REAL> KabsoluteInv(TPZManVector<REAL> & x);
    
    /** @brief Set the material indexes */
    void SetMatIDs(TPZVec<int> &matids) {fmaterialIds=matids;}
    
    /** @brief Get the material indexes */
    TPZVec<int> GetMatIDs() {return fmaterialIds;}
    
    /** @brief Set initial bcs and bcs */
    void SetBC(TPZStack< TPZVec<REAL> > bcini,TPZStack< TPZVec<REAL> > bc){
        
        if (bc.size() ==0  && bcini.size() == 0) {
            std::cout << "There is not boundaries, you give me = " << bc.size() << std::endl;
            DebugStop();
        }
        
        if (bc[0].size() !=4  && bcini[0].size() != 4) {
            std::cout << "The number of parameter must to be equal 4, you give me = " << bc.size() << std::endl;
            DebugStop();
        }
        
        fInitial_BCs = bcini;
        fBCs = bc;
    }
    
    /** @brief Get initial bcs */
    TPZStack< TPZVec<REAL> > GetInitialBC(){
        return fInitial_BCs;
    }
    
    /** @brief Get bcs */
    TPZStack< TPZVec<REAL> > GetBC(){
        return fBCs;
    }
    
};

#endif