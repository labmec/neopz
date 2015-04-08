#include "TPZMatComplexExample.h"
#include "pzbndcond.h"

TPZMatComplexExample::TPZMatComplexExample(int id, REAL l, REAL t, REAL eO, STATE (& ur)( TPZVec<REAL>,REAL),STATE (& er)( TPZVec<REAL>,REAL)) : TPZMaterial(id), fUr(ur), fEr(er), fLambda(l),fTheta(t),fEZero(eO)
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

TPZMatComplexExample::TPZMatComplexExample(int id) : TPZMaterial(id), fUr(urDefault),
fEr(erDefault), fLambda(1e10-9),fTheta(0),fEZero(1)
{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
  
}

/** @brief Default constructor */
TPZMatComplexExample::TPZMatComplexExample() : TPZMaterial(), fUr(urDefault),
fEr(erDefault), fLambda(1e10-9),fTheta(0),fEZero(1)

{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
  
}


TPZMatComplexExample::TPZMatComplexExample(const TPZMatComplexExample &mat) : TPZMaterial(mat), fUr(urDefault),
fEr(erDefault), fLambda(1e10-9),fTheta(0),fEZero(1)
{
  fW=2.*M_PI*M_C/fLambda;
  fKZero=fW*sqrt(M_EZERO*M_UZERO);
    
}

TPZMatComplexExample::~TPZMatComplexExample()
{
    
}
//BATE COM A FORMULAÇÃO FRACA
void TPZMatComplexExample::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL> &x = data.x;
    int nshape=phi.Rows();
  

    for(int i = 0 ; i<nshape ; i++)
    {
        const STATE rhs = 0*phi(i,0);
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
            const STATE stiff=phi(i,0)*phi(j,0)*fKZero*fKZero*(fEr(x,fLambda)-(1./fUr(x,fLambda))*sin(fTheta)*sin(fTheta))-(1./fUr(x,fLambda))*dphi(0,i)*dphi(0,j);
            ek(i,j) += stiff*weight;
        }
    }
    
    
}

void TPZMatComplexExample::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    TPZFMatrix<REAL> &phi = data.phi;
    TPZVec<REAL> &x = data.x;
    int nshape=phi.Rows();
    
    REAL BIG = TPZMaterial::gBigNumber;//sera posto na matriz K
    STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
    STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
    BIG=BIG*BIG;
    switch (bc.Type()) {
        case 0: // Dirichlet
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = phi(i,0)*BIG*v2;
                ef(i,0) += rhs*weight;
                for(int j=0;j<nshape;j++)
                {
                    const STATE stiff = phi(i,0)*phi(j,0);
                    ek(i,j) += stiff*weight*BIG;
                }
            }
            break;
            
        case 1: // Neumann
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = 0*phi(i,0);
                ef(i,0) += rhs*weight+v2;
            }
            break;
            
        case 2: // Mista
            for(int i = 0 ; i<nshape ; i++)
            {
              ef(i,0) += weight+v2*phi(i,0);
                for(int j=0;j<nshape;j++)
                {
                    ek(i,j) += v1*phi(i,0)*phi(j,0);
                }
            }
            break;
            
        default:
            DebugStop();
            break;
    }
    
}


STATE urDefault( TPZVec<REAL> x ,REAL l)
{
  return 1;
}

STATE erDefault( TPZVec<REAL> x, REAL l)
{
  return 1;
}









