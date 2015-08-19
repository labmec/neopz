#include "TPZMatStripline.h"
#include "pzbndcond.h"

TPZMatStripline::TPZMatStripline(int id, REAL freq, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>)) : TPZVecL2(id), fUr(ur), fEr(er), fFreq(freq)
{
  
}

TPZMatStripline::TPZMatStripline(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fFreq(1.0e9)
{
  
}

/** @brief Default constructor */
TPZMatStripline::TPZMatStripline() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fFreq(1.0e9)
{
  
}


TPZMatStripline::TPZMatStripline(const TPZMatStripline &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr), fFreq(mat.fFreq)
{
  
}

TPZMatStripline::~TPZMatStripline()
{
    
}

void TPZMatStripline::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  TPZManVector<STATE,3> force(3);
  if(fForcingFunction) {
    fForcingFunction->Execute(data.x,force);
  }
  
  
  // Setting the phis
  TPZFMatrix<REAL> &phiQ = data.phi;
  
  int phrq;
  phrq = data.fVecShapeIndex.NElements();
  
  
  //Calculate the matrix contribution for flux. Matrix A
  for(int iq=0; iq<phrq; iq++)
  {
    //ef(iq, 0) += 0.;
    int ivecind = data.fVecShapeIndex[iq].first;
    int ishapeind = data.fVecShapeIndex[iq].second;
    TPZFNMatrix<3,REAL> ivec(3,1,0.);
    for(int id=0; id<3; id++){
      ivec(id,0) = data.fNormalVec(id,ivecind);
    }
    STATE ff = 0.;
    for (int i=0; i<3; i++) {
      ff += ivec(i,0)*force[i];
    }
    
    ef(iq,0) += weight*ff*phiQ(ishapeind,0);
    
    for (int jq=0; jq<phrq; jq++)
    {
      TPZFNMatrix<3,REAL> jvec(3,1,0.);
      int jvecind = data.fVecShapeIndex[jq].first;
      int jshapeind = data.fVecShapeIndex[jq].second;
      
      for(int id=0; id<3; id++){
        jvec(id,0) = data.fNormalVec(id,jvecind);
      }
      
      //jvecZ.Print("mat1 = ");
      REAL prod1 = ivec(0,0)*jvec(0,0) + ivec(1,0)*jvec(1,0) + ivec(2,0)*jvec(2,0);
      ek(iq,jq) += weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
    }
  }
}


void TPZMatStripline::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatStripline::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatStripline::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatStripline::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatStripline::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatStripline::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}


STATE urDefault( TPZVec<REAL> x )
{
  return 1.0;
}

STATE erDefault( TPZVec<REAL> x )
{
  return 1.0;
}









