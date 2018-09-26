//$Id: pzelastoplasticmem.h,v 1.7 2009-10-04 05:44:22 erick Exp $

#ifndef PZELASTOPLASTICMEM_H
#define PZELASTOPLASTICMEM_H

#include "TPZMaterial.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"

/**
 * This class defines the material memory for a standar elastoplastic calculation.
 */
class TPZElastoPlasticMem : public TPZSavable
{
    
private:
    
public:
    
    virtual int ClassId() const;
    
    TPZElastoPlasticMem();
    
    TPZElastoPlasticMem(const TPZElastoPlasticMem & other);
    
    const TPZElastoPlasticMem & operator = (const TPZElastoPlasticMem & other);
    
    virtual ~TPZElastoPlasticMem();
    
    const std::string Name() const;
    
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    void Read(TPZStream &buf, void *context);
    
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
    
//    /// Set elastoplastic strain state
//    virtual void SetElastoPlasticState(const TPZPlasticState<REAL> & elastoplastic_state);
//    
//    /// Set elastoplastic strain state
//    virtual TPZPlasticState<REAL> & GetElastoPlasticState();
//    
//    /// Set stress state
//    virtual void SetSigma(const TPZTensor<REAL> & sigma);
//    
//    /// Set stress state
//    virtual TPZTensor<REAL> & GetSigma();
//    
//    /// Set displacement field
//    virtual void SetDisplacement(const TPZManVector<REAL,3> & u);
//    
//    /// Set displacement field
//    virtual TPZManVector<REAL,3> & GetDisplacement();
//    
//    /// Set number of plastic steps
//    virtual void SetPlasticSteps(const int & plastic_steps);
//    
//    /// Set number of plastic steps
//    virtual int & GetPlasticSteps();
//    
//    /// Set yield function value
//    virtual void SetPhi(const REAL & phi);
//    
//    /// Set yield function value
//    virtual REAL & GetPhi();
    
};

#endif

