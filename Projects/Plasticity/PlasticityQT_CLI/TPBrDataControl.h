#ifndef TPBrDataControl_H
#define TPBrDataControl_H

#include <iostream>
#include <fstream>
#include "TPBrLaboratoryData.h"
#include "pzvec.h"

#include "TPZSandlerDimaggio.h"

class TPBrDataControl
{
  
private:
  
    /// modelo plastico que esta na interface grafica
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSandler;
    
    /// ultimo indice que foi utilizado
    int fCounter;
    
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
    
    void DeleteLabData (int med_idx) {
      if (fMedicoes.find(med_idx) == fMedicoes.end())
	DebugStop();
      fMedicoes[med_idx].DeleteAllSimulations();
      fMedicoes.erase(med_idx);
    }
    
    int InsertLaboratoryData (const TPBrLaboratoryData &labdataobj) {
      int labdataidx = GenerateNewIndex();
      fMedicoes[labdataidx] = labdataobj;
      return labdataidx;
    }
    
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
