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

    void DeleteGlobalId(int globalid)
    {
        if(fMapSimMed.find(globalid) != fMapSimMed.end())
        {
            fMapSimMed.erase(globalid);
        }
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

    //    void GetMed (int medid, TPBrLaboratoryData &labdataobj) {
    //      if (fMedicoes.find(medid) == fMedicoes.end()) //nao eh medicao
    //            DebugStop();

    //      labdataobj = fMedicoes[medid];
    //    }

    //    void GetSim (int simid, TPBrSimulationData &simdataobj) {
    //      if (fMapSimMed.find(simid) == fMapSimMed.end()) //nao existe essa simulacao
    //            DebugStop();

    //      int medid = fMapSimMed[simid];

    //      std::cout << "(GetSim method)MedGID: " << fMedicoes[medid].GlobalId() << std::endl;
    //      std::cout << "(GetSim method)SimGID: " << fMedicoes[medid].fSimulacoes[simid].GlobalId() << std::endl;
    //      std::cout << "(GetSim method)fMedicoes[medid]: " << fMedicoes[medid].fSimulacoes.size() << std::endl;
    //      simdataobj = fMedicoes[medid].fSimulacoes[simid];
    //    }

    //    int isMed (int xid) {

    //        std::map<int,TPBrLaboratoryData>::iterator it;

    //        for(it = fMedicoes.begin(); it != fMedicoes.end(); ++it) {
    //            if ((*it).first == xid)
    //                return 1;
    //        }

    //        return 0;
    //    }

    //    int GetMedId (int simid){
    //        if (fMapSimMed.find(simid) == fMapSimMed.end()) //nao existe essa simulacao
    //            DebugStop();

    //        int medid = fMapSimMed[simid];

    //        return medid;
    //    }
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
