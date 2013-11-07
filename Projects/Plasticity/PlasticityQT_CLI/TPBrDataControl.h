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

public:
    TPBrDataControl();
    
    int OpenLabFile(const std::string &filename);
    
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
    
    std::map<int, TPBrLaboratoryData> fMedicoes;
    
    
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
