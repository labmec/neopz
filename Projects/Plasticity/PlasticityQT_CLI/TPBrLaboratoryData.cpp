#include "TPBrLaboratoryData.h"
#include "TPBrPlasticitySimulation.h"

TPBrLaboratoryData::TPBrLaboratoryData()
{
    fstart_idx = -1;
    fend_idx = -1;
}

int TPBrLaboratoryData::RunSimulation () {
  
  // RUN SIMULATION AND RETURN ITS INDEX
  
  TPBrPlasticitySimulation newSimulation;
  newSimulation.SetSimulationInitialStep(fstart_idx);
  //newSimulation.SetSandlerDimaggio();
  
  
  int pos = fSimulacoes.size()+1;
//  fSimulacoes.Resize(pos);
  
  pos = pos - 1;
  //Simulacoes[pos] = newSimulation.GetSimulationData();
  
  //return inserted position
  return pos;
}