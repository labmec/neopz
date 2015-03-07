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
    
    /** @brief Characteristic Pressure - Pa */
    fPref=0;
    
    /** @brief Characteristic Density - kg/m3 */
    fRhoref=0;
    
    /** @brief Characteristic viscosity - Pa s */
    fEtaref=0;
    
    /** @brief Density at P of reference - kg/m3 */
    fRhoref=0;
    
    /** @brief Porosity at P of reference - */
    fPhiref=0;
    
    /** @brief Rock Compressibility 1/pa - */
    fcrock=0;
    
    /** @brief Fluid Compressibility 1/pa - */
    fcfluid=0;
    
    /** @brief Material indexes */
    fmaterialIds.Resize(0);
    
    /** @brief absolute permeability */
    fKab.Resize(2,2);
    fKab.Zero();
    
    /** @brief absolute permeability inverse */
    fKabinv.Resize(2,2);
    fKabinv.Zero();
	
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
    poros = fPhiref*(1.0+fcrock*(P-fPref));
    dPorosDp = fPhiref*fcrock;
}

/**
 * @brief \f$ Oil density RhoOil = RhoOil( P ) \f$
 * @param P fluid pressure
 */
void ReservoirData::Density(REAL P, REAL &rho, REAL &drhoDp)
{
    rho = fRhoref*(1.0+fcfluid*(P-fPref));
    drhoDp = fRhoref*fcfluid;
}

/**
 * @brief Oil viscosity. \f$ OilViscosity = ViscOil( P ) \f$
 * @param P fluid pressure
 */
void ReservoirData::Viscosity(REAL P, REAL &eta, REAL &detaDp)
{
    eta = fEtaref;
    detaDp = 0.0;
}
