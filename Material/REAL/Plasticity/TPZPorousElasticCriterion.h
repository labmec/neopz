//
//  TPZPorousElasticCriterion.h
//  pz
//
//  Created by Omar Durán on 1/16/19.
//

#ifndef TPZPorousElasticCriterion_h
#define TPZPorousElasticCriterion_h

#include <stdio.h>
#include "TPZPlasticState.h"
#include "TPZPlasticStep.h"
#include "TPZPorousElasticResponse.h"

class TPZPorousElasticCriterion : public TPZPlasticBase, public TPZPlasticCriterion
{
    
public:
    
    enum {NYield=3};
    
public:
    
    typedef TPZPorousElasticCriterion m_NYields;
    
    REAL m_kappa = 0.0;
    
    TPZPlasticState<STATE> m_N;
    
    TPZPorousElasticResponse m_ER;
    
public:
    
    TPZPorousElasticCriterion();
    
    TPZPorousElasticCriterion(const TPZPorousElasticCriterion &cp);
    
    TPZPorousElasticCriterion & operator=(const TPZPorousElasticCriterion &cp);
    
    void Read(TPZStream& buf, void* context) override;
    
    void Write(TPZStream& buf, int withclassid) const override;
    
    virtual int IntegrationSteps()const override;
    
    virtual const char * Name() const override
    {
        return "TPZPorousElasticCriterion";
    }
    
    virtual void Print(std::ostream & out) const override
    {
        out << "Classe: " << this->Name();
        m_N.Print(out);
    }
    
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;

    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent = NULL) override;
    
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) override {
        std::cerr << "Deprecated gradient calculation is incorporated on ApplyStrainComputeSigma method." << std::endl;
        DebugStop();
    }
    
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) override;
    
    virtual TPZPlasticState<REAL> GetState() const override;

    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override;
    
    virtual void SetElasticResponse(TPZPorousElasticResponse &ER) override
    {
        m_ER = ER;
    }
    
    virtual TPZPorousElasticResponse GetElasticResponse() const override
    {
        return m_ER;
    }

    virtual void SetState(const TPZPlasticState<REAL> &state) override;
    
    
    virtual int ClassId() const override;
    
    TPZPlasticCriterion& GetYC() override{
        return *this;
    }
    
    void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        TPZTensor<STATE> sigmaTensor;
        sigmaTensor.XX() = sigma[0];
        sigmaTensor.YY() = sigma[1];
        sigmaTensor.ZZ() = sigma[2];
        Phi(sigmaTensor, yield);
    }
    
    virtual int GetNYield() const override {
        return as_integer(NYield);
    }
    
};


#endif /* TPZPorousElasticCriterion_h */
