#include "pzlog.h"
#include "TPZFracData.h"

#include <iostream>
#include <cmath>


TPZFracData::TPZFracData()
{
  fDeltaT = 0.0;
    fRho = 0.0;
    fPhi = 0.0;
  fTime = 0.;
  fTheta = 1.0;
  this->SetCurrentState();
}


TPZFracData::~TPZFracData()
{
  
}

/**
 * @brief \f$ Rock porosity function. \f$ Phi = Phi( p ) \f$
 * @param p pressure
 */
void TPZFracData::Porosity(REAL p, REAL &porosity, REAL &dPorosityDp)
{
    const REAL Rockcomp = (1.0e-8);
    const REAL pref = (10.0e6);
    porosity = fPhi*exp(Rockcomp*(p-pref));
    dPorosityDp = fPhi*Rockcomp*exp(Rockcomp*(p-pref));
}

/**
 * @brief \f$ Fluid density function. \f$ RhoFluid = RhoFluid( p ) \f$
 * @param p pressure
 */
void TPZFracData::Density(REAL p, REAL &RhoFluid, REAL &dRhoDp)
{
    const REAL Oilcomp = (1.0e-8);
    const REAL pref = (10.0e6);
    RhoFluid = fRho*exp(Oilcomp*(p-pref));
    dRhoDp = fRho*Oilcomp*exp(Oilcomp*(p-pref));
}

/**
 * @brief Fluid viscosity function. \f$ FluidViscosity = Visc( p ) \f$
 * @param p pressure
 */
void TPZFracData::Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp)
{
    // Constant for a while.
    FluidViscosity = fmu;
    dFluidViscosityDp = 0;
}

/** @brief Returns absolute permeability. */
TPZFMatrix<STATE> TPZFracData::K(){
    return fKab;
}

/** @brief Returns absolute permeability inverse. */
TPZFMatrix<STATE> TPZFracData::Kinv(){
    TPZFMatrix<STATE> Kinverse(2,2,0.0);
    REAL Constant = (-1.0*fKab(0,1)*fKab(1,0)+fKab(0,0)*fKab(1,1));

    if (fKab.Rows() != 2) {
        std::cout << "Absolute permeability must to be 2 Rows 2 Cols Matrix" << std::endl;
        DebugStop();
    }
        
    Kinverse(0,0) =     1.0*fKab(1,1)/(Constant);
    Kinverse(0,1) = -   1.0*fKab(0,1)/(Constant);
    Kinverse(1,0) = -   1.0*fKab(1,0)/(Constant);
    Kinverse(1,1) =     1.0*fKab(0,0)/(Constant);

    return Kinverse;
}

REAL TPZFracData::G()
{
#ifdef DEBUG
    if (fE < 0. || fnu < 0.) {
        PZError << "ERROR - Elasticity or poisson not set!" << std::endl;
        DebugStop();
    }
#endif
    return fE / (2.*(1+fnu));
}

void TPZFracData::SetNextTime()
{
    if(fTime + fDeltaT > fTtot)
    {
        fDeltaT = fTtot - fTime;
    }
    fTime += fDeltaT;
    std::cout << "\n----------------- Actual Time of Simulation = " << fTime << " -----------------" << std::endl;
}

void TPZFracData::SetPostProcessFileName(std::string &postProcessFileName)
{
    fpostProcessFileName = postProcessFileName;
}

std::string TPZFracData::PostProcessFileName()
{
    return fpostProcessFileName;
}