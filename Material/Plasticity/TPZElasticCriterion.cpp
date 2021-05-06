#include "TPZElasticCriterion.h"

TPZElasticCriterion::TPZElasticCriterion()
{
  
}

TPZElasticCriterion::TPZElasticCriterion(const TPZElasticCriterion &other) : fN(other.fN), fER(other.fER)
{
}

TPZElasticCriterion & TPZElasticCriterion::operator=(const TPZElasticCriterion &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    fN = other.fN;
    fER = other.fER;
    return *this;
}

void TPZElasticCriterion::Read(TPZStream& buf, void* context) {
    fN.Read(buf, context);
    fER.Read(buf, context);
}

void TPZElasticCriterion::Write(TPZStream& buf, int withclassid) const {
    fN.Write(buf, withclassid);
    fER.Write(buf, withclassid);
}

void TPZElasticCriterion::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent)
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
    
    fER.ComputeStress(epsTotal, sigma);
    fN.m_eps_t = epsTotal;
    
    if (require_tangent_Q) {
        const REAL lambda = fER.Lambda();
        const REAL mu = fER.G();
        
        // Linha 0
        tangent->PutVal(_XX_,_XX_, lambda + 2. * mu);
        tangent->PutVal(_XX_,_YY_, lambda);
        tangent->PutVal(_XX_,_ZZ_, lambda);
        
        // Linha 1
        tangent->PutVal(_XY_,_XY_, 2. * mu);
        
        // Linha 2
        tangent->PutVal(_XZ_,_XZ_, 2. * mu);
        
        // Linha 3
        tangent->PutVal(_YY_,_XX_, lambda);
        tangent->PutVal(_YY_,_YY_, lambda + 2. * mu);
        tangent->PutVal(_YY_,_ZZ_, lambda);
        
        // Linha 4
        tangent->PutVal(_YZ_,_YZ_, 2. * mu);
        
        // Linha 5
        tangent->PutVal(_ZZ_,_XX_, lambda);
        tangent->PutVal(_ZZ_,_YY_, lambda);
        tangent->PutVal(_ZZ_,_ZZ_, lambda + 2. * mu);
    }
    
}

void TPZElasticCriterion::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{
  
  std::cout<< " \n this method is not implemented in TPZElasticCriteria. ";
  DebugStop();
  
}

void TPZElasticCriterion::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal)
{
    TPZFNMatrix<36,REAL> Dep(6,6,0.0);
    TPZTensor<REAL> eps(0.),sigma;
    ApplyStrainComputeSigma(eps, sigma, &Dep);
    TPZFNMatrix<6,REAL> stressmat(6,1);
    stressmat(_XX_) = GivenStress[_XX_];
    stressmat(_YY_) = GivenStress[_YY_];
    stressmat(_XY_) = GivenStress[_XY_];
    stressmat(_ZZ_) = GivenStress[_ZZ_];
    stressmat(_XZ_) = GivenStress[_XZ_];
    stressmat(_YZ_) = GivenStress[_YZ_];
    Dep.Solve_LDLt(&stressmat);
    epsTotal[_XX_] = stressmat(_XX_);
    epsTotal[_YY_] = stressmat(_YY_);
    epsTotal[_XY_] = stressmat(_XY_);
    epsTotal[_XZ_] = stressmat(_XZ_);
    epsTotal[_YZ_] = stressmat(_YZ_);
    epsTotal[_ZZ_] = stressmat(_ZZ_);
    fN.m_eps_t = epsTotal;
}


TPZPlasticState<STATE>  TPZElasticCriterion::GetState() const
{
  return fN;
}


void TPZElasticCriterion::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{

  phi.resize(3);
  for (int i = 0; i < 3; i++) {
    phi[i] = 0.;
  }
}

void TPZElasticCriterion::SetState(const TPZPlasticState<REAL> &state)
{
  fN = state;
}

int TPZElasticCriterion::IntegrationSteps() const
{
  return 1;
}

int TPZElasticCriterion::ClassId() const{
    return Hash("TPZElasticCriterion") ^ TPZPlasticBase::ClassId() << 1;
}
