#include "TPBrDataControl.h"
#include <iostream>

#include "pzlog.h"

int main(int argc, char *argv[])
{
  InitializePZLOG();
      string FileName, dirname = PZSOURCEDIR;
      FileName = dirname + "/Projects/Plasticity/";
      FileName += "ensaio_all_columns.txt";

      REAL poisson, E, A, B, C, R, D, W;
      E = 29269;
      poisson = 0.203;
      A = 120; //616.67;
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
      
      DADOS.Set_Med_start_idx(med_idx, 600);//100);
      
      DADOS.Set_Med_elastic_trans_idx(med_idx, 750);
			
      
      TPBrStrainStressDataBase *basedata = DADOS.getObj(med_idx);
      TPBrLaboratoryData *labdata = dynamic_cast<TPBrLaboratoryData *>(basedata);
			labdata->IdentifyElasticity(E,poisson);
			sandlerObj.fER.SetUp(E,poisson);
      if(!labdata) DebugStop();
      
      std::cout << "VAI SIMULAR: Start idx: " << labdata->Get_start_idx() << " End idx: " << labdata->Get_end_idx() << std::endl;

      int idx_sim = labdata->RunSimulation(sandlerObj);
      int idx_med = labdata->GlobalId();
      if(idx_med != med_idx) DebugStop();

      std::cout << "SIMULADO: Med idx: " << idx_med << " Sim idx: " << idx_sim << std::endl;
      
      std::cout << "END" << std::endl;

}
