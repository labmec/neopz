#ifndef TPBrDataControl_H
#define TPBrDataControl_H

#include <iostream>
#include <fstream>
#include "TPBrLaboratoryData.h"
#include "TPBrSimulationData.h"
#include "pzvec.h"

#include "TPZSandlerDimaggio.h"

class TPBrDataControl
{
  
private:
  
    /// modelo plastico que esta na interface grafica
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSandler;
    
    /// ultimo indice que foi utilizado
    int fCounter;
    
    /// mapeamento correspondencia entre simulacao e medicao [SimID, MedID]
    std::map<int, int> fMapSimMed;


public:

    std::map<int, TPBrLaboratoryData> fMedicoes;


    TPBrDataControl();
    
    inline void SetSandlerDimaggio(const TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &copy)
    {
        fSandler = copy;
    }
    
    inline void GetSandlerDimaggio (TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj) {
      obj = fSandler;
    }
    
    int GenerateNewIndex()
    {
        fCounter++;
        return fCounter-1;
    }
    
    void DeleteLabData (int medid) {
      if (fMedicoes.find(medid) == fMedicoes.end()) //nao eh medicao
	DebugStop();
      
      std::map<int,TPBrSimulationData>::iterator ii;
      std::vector<int> SimsToDelete;
      for(ii = fMedicoes[medid].fSimulacoes.begin(); ii != fMedicoes[medid].fSimulacoes.end(); ++ii) {
	int simid = (*ii).first;
	SimsToDelete.push_back(simid);
      }
      for(int i = 0; i < SimsToDelete.size(); i++) {
	std::cout << "Deleting Sim id: " << SimsToDelete[i] << std::endl;
	DeleteSimulationData(SimsToDelete[i]);
      }
      fMedicoes.erase(medid);
      std::cout << "Deleted Med id: " << medid << std::endl;
    }
    
    void DeleteSimulationData (int simid) {
      if (fMapSimMed.find(simid) == fMapSimMed.end()) //nao existe essa simulacao
	DebugStop();
      
      int medid = fMapSimMed[simid];
      fMapSimMed.erase(simid);
      fMedicoes[medid].fSimulacoes.erase(simid);
      std::cout << "Deleted Sim id: " << simid << std::endl;
    }
    
    int InsertLaboratoryData (const TPBrLaboratoryData &labdataobj) {
      int labdataidx = GenerateNewIndex();
      fMedicoes[labdataidx] = labdataobj;
      fMedicoes[labdataidx].SetGlobalId(labdataidx);
      std::cout << "Inserted MED id: " <<  labdataidx << std::endl;
      return labdataidx;
    }
    
    int InsertSimulationData (int medid, const TPBrSimulationData &simdataobj) {
      int simid = fMedicoes[medid].InsertSimulationData(simdataobj);
      fMapSimMed[simid] = medid;
      std::cout << "Inserted Sim id: " << simid  << " on MED " << medid << std::endl;
      return simid;
    }

    int GetMedId (TPBrSimulationData &labdataobj) {
      int MedId;
      
      for (int i=0; i<fMedicoes.size(); i++){
	if (fMedicoes[i] == labdataobj)	
	  MedId = i;
      }
      
      return MedId;
    }
    
    int GetSimId (int medid, TPBrLaboratoryData fMedicoes, TPBrSimulationData &simdataobj) {
      int SimId;
      
      for (int i=0; i<fMedicoes[medid].fSimulacoes.size(); i++){
	if (fMedicoes[medid].fSimulacoes[i] == simdataobj)	
	  SimId = i;
      }
      
      return SimId;
    }
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
