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
  
    /// modelo plastico que gerou esta simulacao
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSandler;

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
    
    TPZVec <TPBrLaboratoryData> Medicoes;
    
    
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
