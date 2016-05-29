//
//  TRMRockProperties.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMRockProperties__
#define __PZ__TRMRockProperties__

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"

class TRMRockProperties{
    
public:
    /** @brief Default constructor */
    TRMRockProperties();

    /** @brief Default constructor */
    TRMRockProperties(const TRMRockProperties &copy)
    {
        DebugStop();
    }
    
    /** @brief Default constructor */
    TRMRockProperties &operator=(const TRMRockProperties &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TRMRockProperties();
    

private:
    
    /**
     * @ingroup Characteristic Parameters in SI
     * @brief Define characteristic parameters for flow in porous media.
     * @since June 09, 2015
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
    
    /** @brief Is GID geometry - */
    bool fIsGIDGeometry;
    
    /** @brief absolute permeability */
    TPZFMatrix<REAL> fKab;
    
    /** @brief absolute permeability inverse */
    TPZFMatrix<REAL> fKabinv;

public:
    
    /** @brief Set the characteristic length - m */
    void SetLref(REAL Lref) {fLref = Lref; }
    
    /** @brief Set the characteristic Pressure - Pa */
    void SetPref(REAL Pref) {fPref = Pref;}
    
    /** @brief Characteristic Pressure - Pa */
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
    
    /** @brief Porosity at P of reference - */
    void SetPhiRef(REAL Phiref) {fPhiref = Phiref;}
    
    /** @brief Porosity at P of reference - */
    REAL PhiRef() {return fPhiref;}
    
    /** @brief Set rock compressibility 1/pa- */
    void SetcRock(REAL cr) {fcrock = cr;}
    
    /** @brief Get rock compressibility 1/pa - */
    REAL CRock() {return fcrock;}
    
    /** @brief Material indexes */
    TPZVec<int> fmaterialIds;
    
    
    /**
     * @brief \f$ Rock porosity. \f$ Phi = Phi( P ) \f$
     * @param P fluid pressure
     * @param poros porosity
     * @param dPorosDp porosity derivative
     */
    void Porosity(REAL P, REAL &poros, REAL &dPorosDp);
    
    /** @brief Set the absolute Permeability - m2 */
    void SetKabsolute(TPZFMatrix<REAL> &Kab)
    {
        fKab = Kab;
        fKabinv = Kab;
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

#endif /* defined(__PZ__TRMRockProperties__) */
