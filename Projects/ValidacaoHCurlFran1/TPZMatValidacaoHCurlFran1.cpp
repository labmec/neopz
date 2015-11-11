#include "TPZMatValidacaoHCurlFran1.h"
#include "pzbndcond.h"

TPZMatValidacaoHCurlFran1::TPZMatValidacaoHCurlFran1(int id, REAL freq, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>)) : TPZVecL2(id), fUr(ur), fEr(er), fFreq(freq)
{
  fW=2.*M_PI*fFreq;
  
}

TPZMatValidacaoHCurlFran1::TPZMatValidacaoHCurlFran1(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fFreq(1.0)
{
  fW=2.*M_PI*fFreq;
}

/** @brief Default constructor */
TPZMatValidacaoHCurlFran1::TPZMatValidacaoHCurlFran1() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fFreq(1.0)
{
  fW=2.*M_PI*fFreq;
}


TPZMatValidacaoHCurlFran1::TPZMatValidacaoHCurlFran1(const TPZMatValidacaoHCurlFran1 &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr), fFreq(mat.fFreq)
{
  fW=2.*M_PI*fFreq;
}

TPZMatValidacaoHCurlFran1::~TPZMatValidacaoHCurlFran1()
{
    
}

void TPZMatValidacaoHCurlFran1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
  TPZManVector<REAL,3> &x = data.x;
  const STATE muR =  fUr(x);
  const STATE epsilonR = fEr(x);
  REAL k0 = 2*M_PI*fFreq*sqrt(M_UZERO*M_EZERO);
  TPZManVector<STATE,3> force(3);
  if(fForcingFunction) {
    fForcingFunction->Execute(x,force);
  }
  
  // Setting the phis
  TPZFMatrix<REAL> &phiQ = data.phi;
  TPZFMatrix<REAL> &dphiQdaxes = data.dphix;//ELEMENTO DEFORMADO
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
  
  TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
  for (int i=0; i<3; i++) {
    ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
    ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
  }
  Cross(ax1, ax2, normal);
  std::cout<<"normal: "<<std::endl;
  std::cout<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<std::endl;
  int phrq;
  phrq = data.fVecShapeIndex.NElements();
  for(int iq=0; iq<phrq; iq++)
  {
    int ivecind = data.fVecShapeIndex[iq].first;
    int ishapeind = data.fVecShapeIndex[iq].second;
    
    TPZManVector<STATE,3> ivecHDiv(3), ivecHCurl(3);
    for(int id=0; id<3; id++){
      ivecHDiv[id] = data.fNormalVec(id,ivecind);
    }
    //ROTATE FOR HCURL
    Cross(normal, ivecHDiv, ivecHCurl);
    STATE ff = 0.;
    for (int i=0; i<3; i++) {
      ff += ivecHCurl[i]*force[i];
    }
    
    ef(iq,0) += weight*ff*phiQ(ishapeind,0);
    TPZManVector<STATE,3> curlI(3), gradPhiI(3);
    for (int i = 0; i<dphiQ.Rows(); i++) {
      gradPhiI[i] = dphiQ(i,ishapeind);//ELEMENTO DEFORMADO
    }
    Cross(gradPhiI, ivecHCurl, curlI);
    
    for (int jq=0; jq<phrq; jq++)
    {
      TPZManVector<STATE,3> jvecHDiv(3), jvecHCurl(3);//JA EM XYZ
      int jvecind = data.fVecShapeIndex[jq].first;
      int jshapeind = data.fVecShapeIndex[jq].second;
      
      for(int id=0; id<3; id++){
        jvecHDiv[id] = data.fNormalVec(id,jvecind);
      }
      //ROTATE FOR HCURL
      Cross(normal, jvecHDiv,jvecHCurl);

      TPZManVector<STATE,3> curlJ(3),  gradPhiJ(3);
      
      for (int i = 0; i<dphiQ.Rows(); i++) {
        gradPhiJ[i] = dphiQ(i,jshapeind);
      }
      Cross(gradPhiJ, jvecHCurl, curlJ);
      //eh necessario forcar curlE dot x = j*k0*sin(theta)*Ez
      
      STATE phiIdotphiJ;
      phiIdotphiJ = ( phiQ(ishapeind,0) * ivecHCurl[0] * phiQ(jshapeind,0) * jvecHCurl[0]);
      phiIdotphiJ += ( phiQ(ishapeind,0) * ivecHCurl[1] * phiQ(jshapeind,0) * jvecHCurl[1]);
      phiIdotphiJ += ( phiQ(ishapeind,0) * ivecHCurl[2] * phiQ(jshapeind,0) * jvecHCurl[2]);
      
      
      STATE stiff = (1./muR) * ( curlJ[0] * curlI[0] + curlJ[1] * curlI[1] + curlJ[2] * curlI[2] );
      stiff += -1. * k0 * k0 * epsilonR * phiIdotphiJ;
      
      ek(iq,jq) += weight * stiff;
    }
  }
}


void TPZMatValidacaoHCurlFran1::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran1::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  
  // Setting the phis
  TPZFMatrix<REAL> &phiQ = data.phi;
  
  int nshape=phiQ.Rows();
  REAL BIG = TPZMaterial::gBigNumber;
  const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
  const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
  switch ( bc.Type() )
  {
    case 0:
      for(int i = 0 ; i<nshape ; i++)
      {
        const STATE rhs = phiQ(i,0) * BIG * v2;
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
          const STATE stiff = phiQ(i,0) * phiQ(j,0) * BIG;
          ek(i,j) += stiff*weight;
        }
      }
      break;
    case 1:
      DebugStop();
      break;
    case 2:
      for(int i = 0 ; i<nshape ; i++)
      {
        const STATE rhs = phiQ(i,0) * v2;
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
          const STATE stiff = phiQ(i,0) *  phiQ(j,0) * v1;
          ek(i,j) += stiff*weight;
        }
      }
      break;
  }
  
}

void TPZMatValidacaoHCurlFran1::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

int TPZMatValidacaoHCurlFran1::VariableIndex(const std::string &name)
{
  if(name == "E") {
      return 2;
  }
  else if ( name == "realE")
  {
    return 3;
  }
  return TPZMaterial::VariableIndex(name);
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatValidacaoHCurlFran1::NSolutionVariables(int var)
{
  int nVar = 0;
  switch (var) {
    case 2:
      nVar = 3;
      break;
    case 3:
      nVar = 3;
      break;
    default:
      nVar = TPZMaterial::NSolutionVariables(var);
      break;
  }
  return nVar;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatValidacaoHCurlFran1::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
  TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
  for (int i=0; i<3; i++) {
    ax1[i] = data.axes(0,i);
    ax2[i] = data.axes(1,i);
  }
  //ROTATE FOR HCURL
  Cross(ax1, ax2, normal);
  std::cout<<"normal sol: "<<std::endl;
  std::cout<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<std::endl;
  
  Solout.Resize(3);
  Cross(normal, data.sol[0], Solout);
  //AQUIFRANAPAGAR
  //Solout=data.sol[0];
  switch (var) {
    case 2://E
    {
#ifdef STATE_COMPLEX
      Solout[0] = std::abs(Solout[0]);
      Solout[1] = std::abs(Solout[1]);
      Solout[2] = std::abs(Solout[2]);
#endif
    }
      break;
    case 3://realE
    {
#ifdef STATE_COMPLEX
      Solout[0] = std::real(Solout[0]);
      Solout[1] = std::real(Solout[1]);
      Solout[2] = std::real(Solout[2]);
//      //AQUIFRANAPAGAR
//      TPZManVector<REAL,3> xAtPoint = data.x;
//      if(xAtPoint[0]+xAtPoint[1]>0)
//      {
//        Solout[0] = 0.;
//        Solout[1] = 0.;
//        Solout[2] = 0.;
//      }
#endif
    }
      break;
      
    default:
    {
    }
      break;
  }
}


STATE urDefault( TPZVec<REAL> x )
{
  return 1.0;
}

STATE erDefault( TPZVec<REAL> x )
{
  return 1.0;
}









