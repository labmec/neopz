#ifndef TPBrLaboratoryData_H
#define TPBrLaboratoryData_H

#include "TPBrStrainStressDataBase.h"
#include "TPBrSimulationData.h"

#include "TPZSandlerDimaggio.h"


class TPBrLaboratoryData : public TPBrStrainStressDataBase
{


public:
  
    TPBrLaboratoryData();  
    TPBrLaboratoryData(const std::string &filename);
    
    TPBrLaboratoryData(const TPBrLaboratoryData &copy) : TPBrStrainStressDataBase(copy),
        fstart_idx(copy.fstart_idx), fend_idx(copy.fend_idx), fSimulacoes(copy.fSimulacoes), felastic_idx(copy.felastic_idx)
    {
        
    }
    
    TPBrLaboratoryData &operator=(const TPBrLaboratoryData &copy)
    {
        TPBrStrainStressDataBase::operator=(copy);
        fstart_idx = copy.fstart_idx;
        fend_idx = copy.fend_idx;
        fSimulacoes = copy.fSimulacoes;
        felastic_idx = copy.felastic_idx;
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
    inline void Set_elastic_trans_idx (int elastic_idx) {
        felastic_idx = elastic_idx;
    }

    virtual int Get_start_idx() const {
        return fstart_idx;
    }
    virtual int Get_end_idx() const {
        return fend_idx;
    }
    virtual int Get_elastic_trans_idx () const {
        return felastic_idx;
    }

    /// Compute the elasticity constants as a function of the stress evolution between fPoreClosureIndex and fElasticTransition
   
    void IdentifyElasticity (REAL &Young, REAL &Poisson);
		
    inline int SizeSimData() const
    {
        return fSimulacoes.size();
    }

    inline int SizeMed() const {
        return this->fSig_Ax.size();
    }

    const TPBrSimulationData *GetSimulation(int globalid) const
    {
        std::map<int, TPBrSimulationData>::const_iterator it;
        it = fSimulacoes.find(globalid);
        if(it == fSimulacoes.end())
        {
            DebugStop();
        }
        else
        {
            return &(it->second);
        }
        return 0;
    }

    TPBrSimulationData *GetSimulation(int globalid)
    {
        std::map<int, TPBrSimulationData>::iterator it;
        it = fSimulacoes.find(globalid);
        if(it == fSimulacoes.end())
        {
            DebugStop();
        }
        else
        {
            return &(it->second);
        }
        return 0;
    }


    void DeleteSimulation(int globalid)
    {
        std::map<int, TPBrSimulationData>::iterator it;
        it = fSimulacoes.find(globalid);
        if(it == fSimulacoes.end())
        {
	  this->Print();
            DebugStop();
        }
        fSimulacoes.erase(globalid);
    }
    
    /// read the input strain and stress from the laboratory file
    void ReadInputStrainStress(const std::string &filename);

    int RunSimulation (TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj);
    
    int InsertSimulationData (const TPBrSimulationData &simdataobj);

    void Print() const
    {
        TPBrStrainStressDataBase::Print();
        std::cout << "fstart_idx " << fstart_idx << std::endl;
        std::cout << "fend_idx " << fend_idx << std::endl;
        std::cout << "felastic_idx " << felastic_idx << std::endl;
        std::map<int, TPBrSimulationData>::const_iterator it;
        for(it = fSimulacoes.begin(); it != fSimulacoes.end(); it++)
        {
            it->second.Print();
        }
    }

protected:
    int fstart_idx;
    int fend_idx;
    int felastic_idx;

    //contains all simulations related to 'this' lab file
    std::map<int, TPBrSimulationData> fSimulacoes;

    
};

#endif // TPBrLaboratoryData_H
