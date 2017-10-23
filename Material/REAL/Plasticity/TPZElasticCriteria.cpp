#include "TPZElasticCriteria.h"

TPZElasticCriteria::TPZElasticCriteria()
{
  
}

TPZElasticCriteria::TPZElasticCriteria(const TPZElasticCriteria &cp) : fN(cp.fN), fER(cp.fER)
{
}

TPZElasticCriteria & TPZElasticCriteria::operator=(const TPZElasticCriteria &cp)
{
    fN = cp.fN;
    fER = cp.fER;
    return *this;
}


void TPZElasticCriteria::Read(TPZStream &buf) { //ok
    fN.Read(buf);
    fER.Read(buf);
}

void TPZElasticCriteria::Write(TPZStream &buf) const { //ok
    fN.Write(buf);
    fER.Write(buf);
}


void TPZElasticCriteria::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{
  fER.Compute(epsTotal, sigma);
}


void TPZElasticCriteria::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
{
  fER.Compute(epsTotal, sigma);
    fN.fEpsT = epsTotal;
  const REAL lambda = fER.Lambda();
  const REAL mu = fER.G();
  Dep.Redim(6,6);
  
  // Linha 0
  Dep(_XX_,_XX_) = lambda + 2. * mu;
  Dep(_XX_,_YY_) = lambda;
  Dep(_XX_,_ZZ_) = lambda;
  
  // Linha 1
  Dep(_XY_,_XY_) = 2. * mu;
  
  // Linha 2
  Dep(_XZ_,_XZ_) = 2. * mu;
  
  // Linha 3
  Dep(_YY_,_XX_) = lambda;
  Dep(_YY_,_YY_) = lambda + 2. * mu;
  Dep(_YY_,_ZZ_) = lambda;

  // Linha 4
  Dep(_YZ_,_YZ_) = 2. * mu;

  // Linha 5
  Dep(_ZZ_,_XX_) = lambda;
  Dep(_ZZ_,_YY_) = lambda;
  Dep(_ZZ_,_ZZ_) = lambda + 2. * mu;
 
}

void TPZElasticCriteria::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{
  
  std::cout<< " \n this method is not implemented in TPZElasticCriteria. ";
  DebugStop();
  
}

void TPZElasticCriteria::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal)
{
    TPZFNMatrix<36,REAL> Dep(6,6);
    TPZTensor<REAL> eps(0.),sigma;
    ApplyStrainComputeDep(eps, sigma, Dep);
    TPZFNMatrix<6,REAL> stressmat(6,1),epsmat(6,1);
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
    fN.fEpsT = epsTotal;
}


TPZPlasticState<STATE>  TPZElasticCriteria::GetState() const
{
  return fN;
}


void TPZElasticCriteria::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{

  phi.resize(3);
  for (int i = 0; i < 3; i++) {
    phi[i] = 0.;
  }
}

void TPZElasticCriteria::SetState(const TPZPlasticState<REAL> &state)
{
  fN = state;
}

int TPZElasticCriteria::IntegrationSteps() const
{
  return 1;
}

int TPZElasticCriteria::ClassId() const{
    return Hash("TPZElasticCriteria") ^ TPZPlasticBase::ClassId() << 1;
}