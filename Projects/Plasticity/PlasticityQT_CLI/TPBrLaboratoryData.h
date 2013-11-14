#ifndef TPBrLaboratoryData_H
#define TPBrLaboratoryData_H

#include "TPBrStrainStressDataBase.h"
#include "TPBrSimulationData.h"

#include "TPZSandlerDimaggio.h"


class TPBrLaboratoryData : public TPBrStrainStressDataBase
{
public:
      //contains all simulations related to 'this' lab file
    std::map<int, TPBrSimulationData> fSimulacoes;
  
    TPBrLaboratoryData();  
    TPBrLaboratoryData(const std::string &filename);
    
    TPBrLaboratoryData(const TPBrLaboratoryData &copy) : TPBrStrainStressDataBase(copy),
        fstart_idx(copy.fstart_idx), fend_idx(copy.fend_idx), fSimulacoes(copy.fSimulacoes)
    {
        
    }
    
    TPBrLaboratoryData &operator=(const TPBrLaboratoryData &copy)
    {
        TPBrStrainStressDataBase::operator=(copy);
        fstart_idx = copy.fstart_idx;
        fend_idx = copy.fend_idx;
        fSimulacoes = copy.fSimulacoes;
        return *this;
    }
    
    virtual ~TPBrLaboratoryData()
    {
        
    }
    
    inline void Set_start_idx(int startidx) {
        fstart_idx = startidx;
    }
    inline void Set_end_idx(int endidx) {
        fend_idx = endidx;
    }
    inline int Get_start_idx() {
        return fstart_idx;
    }
    inline int Get_end_idx() {
        return fend_idx;
    }
    
    /// read the input strain and stress from the laboratory file
    void ReadInputStrainStress(const std::string &filename);

    int RunSimulation (TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj);
    
    int InsertSimulationData (const TPBrSimulationData &simdataobj);

protected:
    int fstart_idx;
    int fend_idx;
    
};

#endif // TPBrLaboratoryData_H
