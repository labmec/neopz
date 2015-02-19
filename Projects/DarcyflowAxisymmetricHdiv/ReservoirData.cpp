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
    fRhoRef=0;
    
    /** @brief Porosity at P of reference - */
    fPhiRef=0;
    
    /** @brief absolute permeability */
    fKab.Resize(2,2);
    fKab.Zero();
	
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
    REAL cR = 1.0e-10;
    poros = fPhiRef*(1.0+cR*(P-fPref));
    dPorosDp = fPhiRef*cR;
}

/**
 * @brief \f$ Oil density RhoOil = RhoOil( P ) \f$
 * @param P fluid pressure
 */
void ReservoirData::RhoOil(REAL P, REAL &Rho, REAL &dRhoDp)
{
    REAL cf = 1.0e-8;
    Rho = fRhoref*(1.0+cf*(P-fPref));
    dRhoDp = fRhoRef*cf;
}

/**
 * @brief Oil viscosity. \f$ OilViscosity = ViscOil( P ) \f$
 * @param P fluid pressure
 */
void ReservoirData::OilViscosity(REAL P, REAL &Viscosity, REAL &dViscosityDp)
{
    Viscosity = 0.005; // 5 [cP] -> 0.005 [Pa s]
    dViscosityDp = 0.0;
}
