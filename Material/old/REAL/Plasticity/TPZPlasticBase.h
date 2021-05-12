/* 
 * File:   TPZPlasticBase.h
 * Author: quinelato
 *
 * Created on 30 de Outubro de 2017, 09:32
 */

#ifndef TPZPLASTICBASE_H
#define TPZPLASTICBASE_H

#include "TPZSavable.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"

class TPZPlasticBase : public TPZSavable {
public:

    int ClassId() const override;

    virtual ~TPZPlasticBase();
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) = 0; //  Candidate to be deprecated.
    // @TODO:: Rename to ComputeStress
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,  TPZFMatrix<REAL> * tangent = NULL) = 0;
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) = 0; //  Candidate to be deprecated.
    // @TODO:: Rename to ComputeStrain
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) = 0;
    
    virtual void SetState(const TPZPlasticState<REAL> &state) = 0;
    virtual TPZPlasticState<REAL> GetState() const = 0;
    
    // @TODO:: Rename to GetYieldCriterion
    virtual TPZPlasticCriterion& GetYC() = 0;
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const = 0;
    virtual int IntegrationSteps() const;
    
    virtual void SetElasticResponse(TPZElasticResponse &ER) = 0;
    virtual TPZElasticResponse GetElasticResponse() const = 0;
    virtual const char * Name()const = 0;
    virtual void Print(std::ostream & out)const = 0;
    virtual void Write(TPZStream& buf, int withclassid) const override = 0;
    virtual void Read(TPZStream& buf, void* context)  override = 0;

    
    virtual void ResetPlasticStrain();
    
};


#endif /* TPZPLASTICBASE_H */

