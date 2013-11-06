#include "TPBrLaboratoryData.h"
#include "../TPZPlasticitySimulation.h"

TPBrLaboratoryData::TPBrLaboratoryData()
{
    fstart_idx = -1;
    fend_idx = -1;
}

int TPBrLaboratoryData::RunSimulation (TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj) {
  
  // RUN SIMULATION AND RETURN ITS INDEX
  
  TPZPlasticitySimulation newSimulation;
  
  newSimulation.ReadInputStrainStress(fSig_Ax, fEps_Ax, fSig_Lat, fEps_Lat);
  newSimulation.SetSimulationInitialStep(fstart_idx);
  newSimulation.SetSandlerDimaggio(obj);
  newSimulation.PerformSimulation();
  
  
  int pos = fSimulacoes.size()+1;
//  fSimulacoes.Resize(pos);
  
  pos = pos - 1;
  //Simulacoes[pos] = newSimulation.GetSimulationData();
  
  //return inserted position
  return pos;
}