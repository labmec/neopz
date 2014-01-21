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

    std::map<int, TPBrLaboratoryData> fMedicoes;

public:



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

    int SizeLabData() const
    {
        return fMedicoes.size();
    }

    int SizeMed (int global_id) {
        return fMedicoes[global_id].SizeMed();
    }
    
    void DeleteLabData (int medid) {
      if (fMedicoes.find(medid) == fMedicoes.end()) //nao eh medicao
        DebugStop();
      
      fMedicoes.erase(medid);
      std::cout << "Deleted Med id: " << medid << std::endl;
    }
    
    
    int InsertLaboratoryData (const TPBrLaboratoryData &labdataobj) {
      int labdataidx = GenerateNewIndex();
      fMedicoes[labdataidx] = labdataobj;
      fMedicoes[labdataidx].SetGlobalId(labdataidx);
      std::cout << "Inserted MED id: " <<  fMedicoes[labdataidx].GlobalId() << std::endl;
      return labdataidx;
    }
    
    int InsertSimulationData (int medid, const TPBrSimulationData &simdataobj) {
      int simid = fMedicoes[medid].InsertSimulationData(simdataobj);
      fMapSimMed[simid] = medid;
      std::cout << "Inserted Sim id: " << simid  << " on MED " << medid << std::endl;
      return simid;
    }

    TPBrStrainStressDataBase *getObj (int global_id)
    {
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) //eh medicao
        {
                    return &(fMedicoes[global_id]);
        }

        if ( (fMapSimMed.find(global_id) != fMapSimMed.end()) ) //eh simulacao
        {
            int medid = fMapSimMed[global_id];

            return fMedicoes[medid].GetSimulation(global_id);
        }
        DebugStop();
        return 0;
    }

    int isMed(int global_id) {
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) //eh medicao
        {
                    return 1;
        }
        if ( (fMapSimMed.find(global_id) != fMapSimMed.end()) ) //eh simulacao
        {
            return 0;
        }
        DebugStop();
        return -1;
    }

    void Set_Med_start_idx (int global_id, int new_idx) {
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) { //existe medicao
            fMedicoes[global_id].Set_start_idx(new_idx);
        }
    }

    void Set_Med_end_idx (int global_id, int new_idx) {
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) { //existe medicao
            fMedicoes[global_id].Set_end_idx(new_idx);
        }
    }

    void Set_Med_elastic_trans_idx (int global_id, int new_idx) {
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) { //existe medicao
            fMedicoes[global_id].Set_elastic_trans_idx(new_idx);
        }
    }

    int Get_Med_start_idx (int global_id){
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) { //existe medicao
            return fMedicoes[global_id].Get_start_idx();
        }
        return -1;
    }

    int Get_Med_end_idx (int global_id){
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) { //existe medicao
            return fMedicoes[global_id].Get_end_idx();
        }
        return -1;
    }

    int Get_Med_elastic_trans_idx (int global_id){
        if ( (fMedicoes.find(global_id) != fMedicoes.end()) ) { //eh medicao
            return fMedicoes[global_id].Get_elastic_trans_idx();
        }
        return -1;
    }

    void DeleteGlobalId(int globalid)
    {
        if(fMapSimMed.find(globalid) != fMapSimMed.end()) //eh simulacao
        {
#ifdef DEBUG
            if ( (fMedicoes.find(globalid) != fMedicoes.end()) )
              DebugStop(); //NAO EXISTE!
#endif
            int med_id = fMapSimMed[globalid];
#ifdef DEBUG    
            if (fMedicoes.find(med_id) == fMedicoes.end())
              DebugStop();
#endif
            fMedicoes[med_id].DeleteSimulation(globalid); //apaga a simulacao
            fMapSimMed.erase(globalid); //remove do mapeamento
            return;
        }
        if ( (fMedicoes.find(globalid) != fMedicoes.end()) ) //eh medicao
        {
                    fMedicoes.erase(globalid); //apaga a medicao (e suas simulacoes)
                    return;
        }
        DebugStop();
    }

    void Print() {
        std::map<int,TPBrLaboratoryData>::iterator it;
        std::cout << "ID DAS MEDICOES:" << std::endl;
        for(it = fMedicoes.begin(); it != fMedicoes.end(); ++it) {
            int idmed = (*it).first;
            std::cout << "MED: " << idmed << std::endl;

            it->second.Print();
        }
    }
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
