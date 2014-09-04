#include "pzlog.h"
#include "TPZFracData.h"

#include <iostream>


TPZFracData::TPZFracData()
{
  fDeltaT = 0.0;
  fTime = 0.;
  fTtot = 0.;
  fmu = -1.;
  fLfrac = -1.;
  fHf = -1.;
  fE = -1.;
  fnu = -1.;
  fPorderPressure = 0;
  fPorderFlow = 1;
  this->SetCurrentState();
}


TPZFracData::~TPZFracData()
{
  
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