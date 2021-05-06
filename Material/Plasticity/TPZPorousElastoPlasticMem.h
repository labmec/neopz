//
//  TPZPorousElastoPlasticMem.h
//  pz
//
//  Created by Omar Durán on 2/11/19.
//

#ifndef TPZPorousElastoPlasticMem_h
#define TPZPorousElastoPlasticMem_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZPorousElasticResponse.h"

/**
 * This class defines the material memory for a porous elastic elastoplastic calculation.
 */
class TPZPorousElastoPlasticMem : public TPZSavable
{
    
private:
    
public:
    
    int ClassId() const override;
    
    TPZPorousElastoPlasticMem();
    
    TPZPorousElastoPlasticMem(const TPZPorousElastoPlasticMem & other);
    
    const TPZPorousElastoPlasticMem & operator = (const TPZPorousElastoPlasticMem & other);
    
    virtual ~TPZPorousElastoPlasticMem();
    
    const std::string Name() const;
    
    void Write(TPZStream &buf, int withclassid) const override;
    
    void Read(TPZStream &buf, void *context) override;
    
    virtual void Print(std::ostream &out = std::cout) const;
    
    friend std::ostream& operator<<( std::ostream& Out, const TPZPorousElastoPlasticMem & s )
    {
        s.Print(Out);
        return Out;
    }
    
    ///  Elastoplastic strain state
    TPZPlasticState<REAL> m_elastoplastic_state;
    
    ///  Stress state
    TPZTensor<REAL> m_sigma;
    
    /// Displacement field
    TPZManVector<REAL,3> m_u;
    
    /// Number of plastic steps
    int m_plastic_steps;
    
    /// Yield function value
    REAL m_phi;
    
    /// Porous Elastoplastic response (It is required when elasti response depends on spatial variables)
    TPZPorousElasticResponse m_ER;
    
};

#endif /* TPZPorousElastoPlasticMem_hpp */
