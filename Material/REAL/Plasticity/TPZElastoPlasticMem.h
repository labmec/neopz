//$Id: pzelastoplasticmem.h,v 1.7 2009-10-04 05:44:22 erick Exp $

#ifndef PZELASTOPLASTICMEM_H
#define PZELASTOPLASTICMEM_H

#include "TPZMaterial.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"

/**
 * This class defines the material memory for a standar elastoplastic calculation.
 */
class TPZElastoPlasticMem : public TPZSavable
{
    
private:
    
public:
    
    int ClassId() const override;
    
    TPZElastoPlasticMem();
    
    TPZElastoPlasticMem(const TPZElastoPlasticMem & other);
    
    const TPZElastoPlasticMem & operator = (const TPZElastoPlasticMem & other);
    
    virtual ~TPZElastoPlasticMem();
    
    const std::string Name() const;
    
    void Write(TPZStream &buf, int withclassid) const override;
    
    void Read(TPZStream &buf, void *context) override;
    
    virtual void Print(std::ostream &out = std::cout) const;
    
    friend std::ostream& operator<<( std::ostream& Out, const TPZElastoPlasticMem & s )
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
    
    /// Elastoplastic response (It is required when elasti response depends on spatial variables)
    TPZElasticResponse m_ER;
    
};

#endif

