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

class TPZPlasticBase : public TPZSavable {
public:

    virtual int ClassId() const;

    virtual ~TPZPlasticBase();
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) = 0;
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZPlasticState<STATE> &plasticState, TPZTensor<REAL> &sigma) = 0;
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZPlasticState<STATE> &plasticState, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) = 0;
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZPlasticState<STATE> &plasticState, TPZTensor<REAL> &epsTotal) = 0;
    virtual TPZPlasticState<STATE> GetInternalState(const TPZPlasticState<REAL> &externalState) const = 0;
    virtual TPZPlasticState<STATE> GetExternalState(const TPZPlasticState<REAL> &internalState) const = 0;
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZPlasticState<STATE> &plasticState, TPZVec<REAL> &phi) const = 0;
    virtual int IntegrationSteps()const;
    virtual void SetElasticResponse(TPZElasticResponse &ER) = 0;
    virtual TPZElasticResponse GetElasticResponse() const = 0;
    virtual const char * Name()const = 0;
    virtual void Print(std::ostream & out)const = 0;
    virtual void Write(TPZStream& buf, int withclassid) const = 0;
    virtual void Read(TPZStream& buf, void* context) = 0;    
};


#endif /* TPZPLASTICBASE_H */

