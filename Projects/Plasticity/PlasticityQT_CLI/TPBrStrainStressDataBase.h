#ifndef TPBrStrainStressDataBase_H
#define TPBrStrainStressDataBase_H

#include "pzreal.h"
#include "pzvec.h"

class TPBrStrainStressDataBase
{

public:
    TPBrStrainStressDataBase();
    TPZVec<REAL> Sig_Ax;
    TPZVec<REAL> Sig_Lat;
    TPZVec<REAL> Eps_Ax;
    TPZVec<REAL> Eps_Lat;

};

#endif // TPBrStrainStressDataBase_H
