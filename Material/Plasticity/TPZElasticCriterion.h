/*
 *  TPZElasticCriterion
 *
 *  Created by Nathan Shauer on 5/4/14.
 *  Copyright 2014 __LabMeC__. All rights reserved.
 *
 */

#ifndef TPZElasticResponse_h
#define TPZElasticResponse_h

#include "TPZPlasticState.h"
#include "TPZPlasticStep.h"
#include "TPZElasticResponse.h"

class TPZElasticCriterion : public TPZPlasticBase, public TPZPlasticCriterion
{
    
public:
    
    enum {NYield=3};
    
    typedef TPZElasticCriterion fNYields;
    
    /// Plastic state
    TPZPlasticState<STATE> fN;
    
    /// Elastic response
    TPZElasticResponse fER;
    
public:
    
    /**
     A unique class identifier
     
     @return The class identifier as integer
     */
    virtual int ClassId() const override;
    
    /**
     Default constructor
     */
    TPZElasticCriterion();
    
    /**
     Copy constructor
     */
    TPZElasticCriterion(const TPZElasticCriterion &cp);
    
    /**
     Assignment constructor
     */
    TPZElasticCriterion & operator=(const TPZElasticCriterion &cp);
    
    
    /**
     Read (persistency)
     
     @param buf a TPZStream object
     @param context pointer to the associated object
     */
    void Read(TPZStream& buf, void* context) override;
    
    
    /**
     Write (persistency)
     
     @param buf a TPZStream object
     @param withclassid the class identifier
     */
    void Write(TPZStream& buf, int withclassid) const override;
    
    /** @brief Return  */
    
    
    /**
     Numer of integration steps

     @return the number of plastic steps in the last load step. Zero indicates elastic loading.
     */
    virtual int IntegrationSteps()const override;
    
    
    /**
     Class name

     @return constant char with the class name
     */
    virtual const char * Name() const override
    {
        return "TPZElasticCriterion";
    }
    
    
    /**
     Print

     @param out ostream object to write the output
     */
    virtual void Print(std::ostream & out) const override
    {
        out << "Classe: " << this->Name();
        fN.Print(out);
    }
    
    /**
     Imposes the specified strain tensor, evaluating the plastic integration if necessary

     @param epsTotal Imposed total strain tensor
     */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;
    
    /**
     Imposes the specified strain tensor and returns the correspondent stress state.

     @param epsTotal Imposed total strain tensor
     @param sigma Stress tensor
     @param De Incremental constitutive relation in Voigt notation (optional)
     */
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * De = NULL) override;
    
    
    /**
     Imposes the specified strain tensor and returns the corresping stress state and tangent matrix. (Deprecated)

     @param epsTotal Imposed total strain tensor
     @param sigma Stress tensor
     @param De Incremental constitutive relation in Voigt notation
     */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &De) override {
        std::cerr << "Deprecated gradient calculation is incorporated on ApplyStrainComputeSigma method." << std::endl;
        DebugStop();
    }
    
    /**
     Computes Strain tensor based on inverse of incremental constitutive relation with a given Stress tensor

     @param sigma Stress tensor
     @param epsTotal Strain tensor
     */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) override;
    

    /**
     Set the plastic state member
     
     @param state The plastic state member
     */
    virtual void SetState(const TPZPlasticState<REAL> &state) override;
    
    /**
     Access to the plastic state member

     @return the plastic state member
     */
    virtual TPZPlasticState<REAL> GetState() const override;
    
    
    /**
     Computes the value of the yield functions for the given strain

     @param epsTotal Strain tensor
     @param phi The valued yield functions
     */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override;
    
    
    /**
     Set the linear elastic response memeber

     @param ER The linear elastic response
     */
    virtual void SetElasticResponse(TPZElasticResponse &ER) override
    {
        fER = ER;
    }
    
    
    /**
     Access to the linear elastic response memeber

     @return The linear elastic response
     */
    virtual TPZElasticResponse GetElasticResponse() const override
    {
        return fER;
    }
    
    
    /**
     Access to the plastic criterion being used

     @return The plastic criterion
     */
    TPZPlasticCriterion & GetYC() override{
        return *this;
    }
    
    
    /**
     Evaluates the yield functions base on a given Stress and hardening

     @param sigma The Stress tensor
     @param hardening The Hardening variable
     @param yield The valued yield functions
     */
    void YieldFunction(const TPZVec<STATE>& sigma, STATE hardening, TPZVec<STATE>& yield) const override{
        TPZTensor<STATE> sigmaTensor;
        sigmaTensor.XX() = sigma[0];
        sigmaTensor.YY() = sigma[1];
        sigmaTensor.ZZ() = sigma[2];
        Phi(sigmaTensor, yield);
    }
    
    /**
     Number of yield surfaces

     @return the number of yield surfaces
     */
    virtual int GetNYield() const override {
        return as_integer(NYield);
    }
    
    
};


#endif /* TPZElasticCriterion_h */
