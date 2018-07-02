/*
 *  TPZElasticCriterion
 *
 *  Created by Nathan Shauer on 5/4/14.
 *  Copyright 2014 __LabMeC__. All rights reserved.
 *
 */

#ifndef TPZELASTICCRITERION
#define TPZELASTICCRITERION

#include "TPZPlasticState.h"
#include "TPZPlasticStep.h"
#include "TPZElasticResponse.h"

class TPZElasticCriterion : public TPZPlasticBase, public TPZPlasticCriterion
{

public:
  enum {NYield=3};
  
public:
	typedef TPZElasticCriterion fNYields;
  
  /** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
  TPZPlasticState<STATE> fN;
  
  TPZElasticResponse fER;
  
public:
  
  /**
   * @brief empty constructor
   */
  TPZElasticCriterion();
  
  /**
   * @brief Copy Constructor
   */
  TPZElasticCriterion(const TPZElasticCriterion &cp);
  
  /**
   * @brief Operator =
   */
  TPZElasticCriterion & operator=(const TPZElasticCriterion &cp);
  
  void Read(TPZStream& buf, void* context) override;
  
  void Write(TPZStream& buf, int withclassid) const override;
  
  /** @brief Return the number of plastic steps in the last load step. Zero indicates elastic loading. */
  virtual int IntegrationSteps()const override;
  
  
  /**
   * @brief Name of the class ina string
   */
  virtual const char * Name() const override
  {
    return "TPZElasticCriteria";
  }
  
  virtual void Print(std::ostream & out) const override
  {
    out << "Classe: " << this->Name();
    fN.Print(out);
  }
    
  /**
   * Imposes the specified strain tensor, evaluating the plastic integration if necessary.
   *
   * @param[in] epsTotal Imposed total strain tensor
   */
  virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;
  
  /**
   * Imposes the specified strain tensor and returns the correspondent stress state.
   *
   * @param[in] epsTotal Imposed total strain tensor
   * @param[out] sigma Resultant stress
   */
  virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent = NULL) override;
  
  
  
  //virtual void ApplySigmaComputeStrain(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal);
  
  /**
   * Imposes the specified strain tensor and returns the corresp. stress state and tangent
   * stiffness matrix.
   *
   * @param[in] epsTotal Imposed total strain tensor
   * @param[out] sigma Resultant stress
   * @param[out] Dep Incremental constitutive relation
   */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) override {
        std::cerr << "Deprecated gradient calculation is incorporated on ApplyStrainComputeSigma method." << std::endl;
        DebugStop();
    }
  
  /**
   * Attempts to compute an epsTotal value in order to reach an imposed stress state sigma.
   * This method should be used only for test purposes because it isn't fully robust. Some
   * materials, specially those perfectly plastic and with softening, may fail when applying
   * the Newton Method on ProcessLoad.
   *
   * @param[in] sigma stress tensor
   * @param[out] epsTotal deformation tensor
   */
  virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) override;
  
  virtual TPZPlasticState<REAL> GetState() const override;
  /**
   * @brief Return the value of the yield functions for the given strain
   * @param[in] epsTotal strain tensor (total strain)
   * @param[out] phi vector of yield functions
   */
  virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override;
  
  virtual void SetElasticResponse(TPZElasticResponse &ER) override
  {
      fER = ER;
  }
  
    virtual TPZElasticResponse GetElasticResponse() const override
    {
        return fER;
    }
  /**
   * @brief Update the damage values
   * @param[in] state Plastic state proposed
   */
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


#endif //TPZELASTICCRITERION
