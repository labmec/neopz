#include "TPBrLaboratoryData.h"
#include "../TPZPlasticitySimulation.h"
#include "TPBrDataControl.h"

TPBrLaboratoryData::TPBrLaboratoryData()
{
    fstart_idx = -1;
    fend_idx = -1;
}

TPBrLaboratoryData::TPBrLaboratoryData(const std::string &filename)
{
    fstart_idx = -1;
    fend_idx = -1;
    ReadInputStrainStress(filename);
}

int TPBrLaboratoryData::InsertSimulationData (const TPBrSimulationData &simdataobj) {
      int simid = DADOS.GenerateNewIndex();
      fSimulacoes[simid] = simdataobj;
      fSimulacoes[simid].SetGlobalId(simid);
      return simid;
}

int TPBrLaboratoryData::RunSimulation (TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj) {
  
  // RUN SIMULATION AND RETURN ITS INDEX
  
  TPZPlasticitySimulation newSimulation;
  
  newSimulation.ReadInputStrainStress(fSig_Ax, fEps_Ax, fSig_Lat, fEps_Lat);
  newSimulation.SetSimulationInitialStep(fstart_idx);
  newSimulation.SetSandlerDimaggio(obj);
  newSimulation.PerformSimulation();
    TPZVec<REAL> sigax,sigr,epsax,epsr;
    newSimulation.GetSimulatedStrainStress(sigax, epsax, sigr, epsr);
  
    
    TPBrSimulationData result;
    int medid = GlobalId();
    result.Set_medicao_idx(medid);
    result.Set_start_idx(fstart_idx);
    result.Set_end_idx(fend_idx);
    result.SetStrainStress(sigax, epsax, sigr, epsr);
    int resultid = DADOS.InsertSimulationData(medid, result);
    
  //return inserted position
  return resultid;
}

/// read the input strain and stress from the laboratory file
void TPBrLaboratoryData::ReadInputStrainStress(const std::string &filename)
{
    std::ifstream input(filename.c_str());
    if (!input) {
        DebugStop();
    }
    int numlines = 0;
    char buf[1024];
    input.getline(buf , 1024);
    STATE x, sig_ax_t, tempo, sig_ax_dev, sig_r, eps_ax, eps_r, eps_v;
    while (input) {
        input >> x >> sig_ax_t >> tempo >> sig_ax_dev >> tempo >> sig_r >> tempo >> eps_ax >> tempo >> eps_r >> tempo >> eps_v;
        if (!input) {
            break;
        }
        if(numlines >= fSig_Ax.size())
        {
            fSig_Ax.Resize(numlines+100, 2);
            fSig_Lat.Resize(numlines+100, 2);
            fEps_Ax.Resize(numlines+100, 2);
            fEps_Lat.Resize(numlines+100, 2);
        }
        fSig_Ax[numlines] = -sig_ax_t;
        fSig_Lat[numlines] = -sig_r;
        fEps_Ax[numlines] = -eps_ax/100.;
        fEps_Lat[numlines] = -eps_r/100.;
        
        //    cout << "i= " << numlines << " " <<  fStressRZInput(numlines,1) << " " <<  fStressRZInput(numlines,0)<< " " << fStrainRZInput(numlines,1)<< " "<< fStrainRZInput(numlines,0) << endl;
        
        numlines++;
    }
    fSig_Ax.Resize(numlines);
    fSig_Lat.Resize(numlines);
    fEps_Lat.resize(numlines);
    fEps_Ax.resize(numlines);
}

