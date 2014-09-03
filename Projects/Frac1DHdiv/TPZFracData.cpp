#include "pzlog.h"
#include "TPZFracData.h"

#include <iostream>


TPZFracData::TPZFracData()
{
  fDeltaT = 0.0;
  fmu = 0.;
  fTime = 0.;
  this->SetCurrentState();
}


TPZFracData::~TPZFracData()
{
  
}

