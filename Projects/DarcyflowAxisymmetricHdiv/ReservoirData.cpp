/*
 *  ReservoirData.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "ReservoirData.h"

ReservoirData::ReservoirData()
{
    /** @brief Characteristic length - m */
    fLref=0;
    
    /** @brief Characteristic Permeability - m2 */
    fKref=0;
    
    /** @brief Porosity at P of reference - */
    fPhiref=0;
    
    /** @brief Characteristic Pressure - Pa */
    fPref=0;
    
    /** @brief Rock Compressibility 1/pa - */
    fcrock=0;
    
    /** @brief Layer average thickness - m */
    fh = 0.0;
    
    /** @brief Layer average radius - m */
    fr = 0.0;
    
    /** @brief well radius - m */
    frw = 0.0;
    
    /** @brief Layer Top depth  - m */
    fDepth = 0.0;
    
    /** @brief Irreducible Saturation of the wetting phase - */
    fS_wett_r = 0.0;
    
    /** @brief Irreducible Saturation of the no wetting phase - */
    fS_nwett_r = 0.0;
    
    
    /** @brief Is GID geometry - */
    fIsGIDGeometry = false;
    
    /** @brief Material indexes */
    fmaterialIds.Resize(0);
    
    /** @brief absolute permeability */
    fKab.Resize(2,2);
    fKab.Zero();
    
    /** @brief absolute permeability inverse */
    fKabinv.Resize(2,2);
    fKabinv.Zero();
    
    /** @brief Absolute permeability vector */
    fKabVec.Resize(0);
    
    /** @brief Definition of the initial boundaries */
    fInitial_BCs.Resize(0);
    
    /** @brief Definition of the boundaries */
    fBCs.Resize(0);
	
}

ReservoirData::~ReservoirData()
{
	
}


/**
 * @brief \f$ Rock porosity. \f$ Phi = Phi( P ) \f$
 * @param P fluid pressure
 */
void ReservoirData::Porosity(REAL P, REAL &poros, REAL &dPorosDp)
{
    poros = fPhiref*(1.0+(fcrock*P-fcrock*fPref));
    dPorosDp = fPhiref*fcrock;
}

/** @brief Get the absolute Permeability - m2 */
TPZFMatrix<REAL> Kabsolute(TPZManVector<REAL> & x){
    
    TPZFMatrix<REAL> k;
    
    return k;
}

/** @brief Get the absolute Permeability inverse - 1/m2 */
TPZFMatrix<REAL> KabsoluteInv(TPZManVector<REAL> & x){
    
    TPZFMatrix<REAL> kinv;
    
    return kinv;
    
}




