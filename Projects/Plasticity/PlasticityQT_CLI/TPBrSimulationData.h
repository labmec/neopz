#ifndef TPBrSimulationData_H
#define TPBrSimulationData_H

#include "TPBrStrainStressDataBase.h"
#include "TPZSandlerDimaggio.h"

class TPBrSimulationData : public TPBrStrainStressDataBase
{

public:
    TPBrSimulationData();
    TPBrSimulationData(int startidx, int endidx, int medidx);

    inline void Set_start_idx(int startidx) {
        start_idx = startidx;
    }
    inline void Set_end_idx(int endidx) {
        end_idx = endidx;
    }
    inline int Get_start_idx(){
        return start_idx;
    }
    inline int Get_end_idx(){
        return end_idx;
    }
    inline void Set_med_idx(int medidx){
        med_idx = medidx;
    }
    inline int Get_med_idx(){
        return med_idx;
    }
    inline void SetUpSandlerDimaggio(REAL poisson, REAL E, REAL A, REAL B, REAL C, REAL R, REAL D, REAL W)
    {
        fSandler.SetUp(poisson, E, A, B, C, R, D, W);
    }


protected:
    int start_idx;
    int end_idx;
    int med_idx;
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSandler;

};

#endif // TPBrSimulationData_H
