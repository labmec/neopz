#ifndef TPBrDataControl_H
#define TPBrDataControl_H

#include <iostream>
#include <fstream>
#include "TPBrLaboratoryData.h"
#include "pzvec.h"

class TPBrDataControl
{

public:
    TPBrDataControl();
    int OpenLabFile(const std::string &filename);
    TPZVec <TPBrLaboratoryData> Medicoes;
};

extern TPBrDataControl DADOS;

#endif // TPBrDataControl_H
