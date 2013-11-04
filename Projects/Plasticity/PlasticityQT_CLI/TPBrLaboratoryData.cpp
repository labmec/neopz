#include "TPBrLaboratoryData.h"
#include "TPBrPlasticitySimulation.h"

TPBrLaboratoryData::TPBrLaboratoryData()
{
    start_idx = -1;
    end_idx = -1;
}

int TPBrLaboratoryData::RunSimulation () {
  
  // RUN SIMULATION AND RETURN ITS INDEX
  
  TPBrPlasticitySimulation newSimulation;
  newSimulation.SetSimulationInitialStep(start_idx);
  //newSimulation.SetSandlerDimaggio();
  
  
  int pos = Simulacoes.size()+1;
  Simulacoes.Resize(pos);
  
  pos = pos - 1;
  //Simulacoes[pos] = newSimulation.GetSimulationData();
  
  //return inserted position
  return pos;
}