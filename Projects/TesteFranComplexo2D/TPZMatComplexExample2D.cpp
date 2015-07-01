#include "TPZMatComplexExample2D.h"
#include "pzbndcond.h"

TPZMatComplexExample2D::TPZMatComplexExample2D(int id, REAL l, REAL t, REAL eO, STATE (& ur)( TPZVec<REAL>,REAL),STATE (& er)( TPZVec<REAL>,REAL)) : TPZVecL2(id), fUr(ur), fEr(er), fLambda(l),fTheta(t),fEZero(eO)
{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
//  std::cout<<"fLambda = "<<fLambda<<"\n";
//  std::cout<<"e0 = "<<M_EZERO<<"\n";
//  std::cout<<"u0 = "<<M_UZERO<<"\n";
//  std::cout<<"c = "<<M_C<<"\n";
//  std::cout<<"fW = "<<fW<<"\n";
//  std::cout<<"fKZero = "<<fKZero<<"\n";
//  std::cout<<" "<<"\n";
  
}

TPZMatComplexExample2D::TPZMatComplexExample2D(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fLambda(1e10-9),fTheta(0),fEZero(1)
{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
  
}

/** @brief Default constructor */
TPZMatComplexExample2D::TPZMatComplexExample2D() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fLambda(1e10-9),fTheta(0),fEZero(1)

{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
  
}


TPZMatComplexExample2D::TPZMatComplexExample2D(const TPZMatComplexExample2D &mat) : TPZVecL2(mat), fUr(urDefault),
fEr(erDefault), fLambda(1e10-9),fTheta(0),fEZero(1)
{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
    
}

TPZMatComplexExample2D::~TPZMatComplexExample2D()
{
    
}

void TPZMatComplexExample2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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


void TPZMatComplexExample2D::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatComplexExample2D::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatComplexExample2D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatComplexExample2D::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatComplexExample2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatComplexExample2D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}


STATE urDefault( TPZVec<REAL> x ,REAL l)
{
  return 1;
}

STATE erDefault( TPZVec<REAL> x, REAL l)
{
  return 1;
}









