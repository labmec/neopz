#include "TPBrDataControl.h"
#include <iostream>

int main(int argc, char *argv[])
{
  
      string FileName, dirname = PZSOURCEDIR;
      FileName = dirname + "/Projects/Plasticity/";
      FileName += "ensaio_all_columns.txt";

      REAL poisson, E, A, B, C, R, D, W;
      E = 29269;
      poisson = 0.203;
      A = 616.67;
      B = 0.0036895;
      C = 111.48;
      D = 0.018768;
      R = 0.91969;
      W = 0.006605;
      
      TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> sandlerObj;
      sandlerObj.SetUp(poisson, E, A, B, C, R, D, W);
      DADOS.SetSandlerDimaggio(sandlerObj);
      
      TPBrLaboratoryData newLabFile (FileName);
      int med_idx = DADOS.InsertLaboratoryData(newLabFile);
      //DADOS.fMedicoes[med_idx].Set_start_idx(100);
      //int idx_sim = DADOS.fMedicoes[med_idx].RunSimulation(sandlerObj);
      
      TPBrLaboratoryData newLabFile2 (FileName);
      int med_idx2 = DADOS.InsertLaboratoryData(newLabFile2);
      //DADOS.fMedicoes[med_idx2].Set_start_idx(100);
      //int idx_sim2 = DADOS.fMedicoes[med_idx2].RunSimulation(sandlerObj);
      //int idx_sim3 = DADOS.fMedicoes[med_idx2].RunSimulation(sandlerObj);
      
      //std::cout << "Simid1: " << idx_sim << " Simid2 " << idx_sim2 << " Simid3 " << idx_sim3  << " Medid1 " << med_idx << " Medid2 " << med_idx2 << std::endl;
      
      //DADOS.DeleteLabData(2);
      //DADOS.DeleteSimulationData(1);
      //DADOS.DeleteLabData(0);
      //DADOS.DeleteLabData(0);
      
      TPBrSimulationData retsim;
      retsim.SetGlobalId(99);
      //DADOS.GetSim(3, retsim);
      
      //std::cout << "---> GetSim " << retsim.GlobalId() << std::endl;
      
      TPBrLaboratoryData retmed;
      //DADOS.GetMed(2, retmed);
      
      //std::cout << "GetMed " << retmed.GlobalId() << " size " << retmed.fSimulacoes.size() << std::endl;
      
      
      std::cout << "END" << std::endl;

}
