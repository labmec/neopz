//
//  TPZPorousElasticCriterion.h
//  pz
//
//  Created by Omar Durán on 1/16/19.
//

#ifndef TPZPorousElasticCriterion_h
#define TPZPorousElasticCriterion_h

#include "TPZPlasticState.h"
#include "TPZPlasticStep.h"
#include "TPZElasticResponse.h"
#include "TPZPorousElasticResponse.h"

class TPZPorousElasticCriterion : public TPZPlasticBase, public TPZPlasticCriterion
{
    
public:
    
    enum {NYield=3};
    
    typedef TPZPorousElasticCriterion fNYields;
    
    /// Elastoplastic state
    TPZPlasticState<STATE> fN;
    
    /// Linearized elastic response
    TPZElasticResponse fER;
    
    /// Porous elastic response
    TPZPorousElasticResponse fPER;
    
public:
    
    /// Default constructor
    TPZPorousElasticCriterion();
    
    /// Copy constructor
    TPZPorousElasticCriterion(const TPZPorousElasticCriterion &cp);
    
    /// Assignment constructor
    TPZPorousElasticCriterion & operator=(const TPZPorousElasticCriterion &cp);
    
    /// Read members
    void Read(TPZStream& buf, void* context) override;
    
    /// Write members
    void Write(TPZStream& buf, int withclassid) const override;

    virtual int IntegrationSteps()const override;
    
    virtual const char * Name() const override
    {
        return "TPZPorousElasticCriterion";
    }
    
    virtual void Print(std::ostream & out) const override
    {
        out << "Classe: " << this->Name();
        fN.Print(out);
        fER.Print(out);
        fPER.Print(out);
    }
    
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;
    
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &eps, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent = NULL) override;
    
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &eps, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) override {
        std::cerr << "Deprecated gradient calculation is incorporated on ApplyStrainComputeSigma method." << std::endl;
        DebugStop();
    }
    
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &eps) override;
    
    virtual TPZPlasticState<REAL> GetState() const override;

    virtual void Phi(const TPZTensor<REAL> &eps, TPZVec<REAL> &phi) const override;
    
    virtual void SetElasticResponse(TPZElasticResponse &ER) override
    {
        fER = ER;
    }
    
    virtual TPZElasticResponse GetElasticResponse() const override
    {
        return fER;
    }
    
    virtual void SetPorousElasticResponse(TPZPorousElasticResponse &PER)
    {
        fPER = PER;
    }
    
    virtual TPZPorousElasticResponse GetPorousElasticResponse() const
    {
        return fPER;
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
