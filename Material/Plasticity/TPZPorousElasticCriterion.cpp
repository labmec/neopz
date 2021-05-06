//
//  TPZPorousElasticCriterion.cpp
//  pz
//
//  Created by Omar Durán on 1/16/19.
//

#include "TPZPorousElasticCriterion.h"

TPZPorousElasticCriterion::TPZPorousElasticCriterion()
{
    
}

TPZPorousElasticCriterion::TPZPorousElasticCriterion(const TPZPorousElasticCriterion &other) : fN(other.fN), fER(other.fER), fPER(other.fPER)
{
}

TPZPorousElasticCriterion & TPZPorousElasticCriterion::operator=(const TPZPorousElasticCriterion &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    fN = other.fN;
    fER = other.fER;
    fPER = other.fPER;
    return *this;
}

void TPZPorousElasticCriterion::Read(TPZStream& buf, void* context) {
    fN.Read(buf, context);
    fER.Read(buf, context);
    fPER.Read(buf, context);
}

void TPZPorousElasticCriterion::Write(TPZStream& buf, int withclassid) const {
    fN.Write(buf, withclassid);
    fER.Write(buf, withclassid);
    fPER.Write(buf, withclassid);
}

//#define LE_Q

void TPZPorousElasticCriterion::ApplyStrainComputeSigma(const TPZTensor<REAL> &eps, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent)
{
    
    bool require_tangent_Q = (tangent != NULL);
    
#ifdef PZDEBUG
    if (require_tangent_Q){
        // Check for required dimensions of tangent
        if (tangent->Rows() != 6 || tangent->Cols() != 6) {
            std::cerr << "Unable to compute the tangent operator. Required tangent array dimensions are 6x6." << std::endl;
            DebugStop();
        }
    }
#endif
    
#ifdef LE_Q
    
    fER.ComputeStress(eps, sigma);
    fN.m_eps_t = eps;

    if (require_tangent_Q) {
        fER.De(*tangent);
    }

    
#else
    
    fPER.ComputeStress(eps, sigma);
    fN.m_eps_t = eps;
    
    if (require_tangent_Q) {
        fPER.De(eps, *tangent);
    }
    
#endif
    
}

void TPZPorousElasticCriterion::ApplyStrain(const TPZTensor<REAL> &eps)
{
    
    std::cout<< " \n this method is not implemented in TPZElasticCriteria. ";
    DebugStop();
    
}

void TPZPorousElasticCriterion::ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsilon)
{
    
#ifdef LE_Q
    
    fER.ComputeStrain(sigma, epsilon);
    fN.m_eps_t = epsilon;
    
    /// Update the linear elastic response
    fER = fPER.EvaluateElasticResponse(epsilon);
    
#else
    
    fPER.ComputeStrain(sigma, epsilon);
    fN.m_eps_t = epsilon;
    
    /// Update the linear elastic response
    fER = fPER.EvaluateElasticResponse(epsilon);
    
#endif
}


TPZPlasticState<STATE>  TPZPorousElasticCriterion::GetState() const
{
    return fN;
}


void TPZPorousElasticCriterion::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{
    phi.resize(3);
    for (int i = 0; i < 3; i++) {
        phi[i] = 0.;
    }
}

void TPZPorousElasticCriterion::SetState(const TPZPlasticState<REAL> &state)
{
    fN = state;
}

int TPZPorousElasticCriterion::IntegrationSteps() const
{
    return 1;
}

int TPZPorousElasticCriterion::ClassId() const {
    return Hash("TPZPorousElasticCriterion") ^ TPZPlasticBase::ClassId() << 1;
}


