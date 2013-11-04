#ifndef TPBrLaboratoryData_H
#define TPBrLaboratoryData_H

#include "TPBrStrainStressDataBase.h"
#include "TPBrSimulationData.h"

class TPBrLaboratoryData : public TPBrStrainStressDataBase
{

public:
    TPBrLaboratoryData();
    
    inline void Set_start_idx(int startidx) {
        start_idx = startidx;
    }
    inline void Set_end_idx(int endidx) {
        end_idx = endidx;
    }
    inline int Get_start_idx() {
        return start_idx;
    }
    inline int Get_end_idx() {
        return end_idx;
    }
    inline void Attribute_sim() {

    }
    int RunSimulation ();

protected:
    int start_idx;
    int end_idx;
    
    //contains all simulations related to 'this' lab file
    TPZVec <TPBrSimulationData> Simulacoes;

};

#endif // TPBrLaboratoryData_H
