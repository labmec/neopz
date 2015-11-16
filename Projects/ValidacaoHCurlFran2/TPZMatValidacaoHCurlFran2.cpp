#include "TPZMatValidacaoHCurlFran2.h"
#include "pzbndcond.h"

TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2(int id, REAL freq, STATE (ur)( TPZVec<REAL>&),STATE (er)( TPZVec<REAL>&), REAL t, REAL scale) : TPZVecL2(id), fUr(ur), fEr(er), fFreq(freq), fScale(scale)
{
  fW=2.*M_PI*fFreq;
  fTheta = t;
  
}

TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fFreq(1.0)
{
  fW=2.*M_PI*fFreq;
  fTheta = 0.;
}

/** @brief Default constructor */
TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fFreq(1.0)
{
  fW=2.*M_PI*fFreq;
  fTheta = 0.;
}


TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2(const TPZMatValidacaoHCurlFran2 &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr), fFreq(mat.fFreq)
{
  fW=2.*M_PI*fFreq;
  fTheta = mat.fTheta;
}

TPZMatValidacaoHCurlFran2::~TPZMatValidacaoHCurlFran2()
{
    
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  TPZManVector<REAL,3> x = data.x;
    for (int i=0; i<3; i++) {
        x[i] /= fScale;
    }
  const STATE muR =  fUr(x);
  const STATE epsilonR = fEr(x);
  REAL k0 = 2*M_PI*fFreq*sqrt(M_UZERO*M_EZERO);
  TPZManVector<STATE,3> force(3);
  if(fForcingFunction) {
    fForcingFunction->Execute(x,force);
  }
  //AQUIFRAN
  // Setting the phis
  TPZFMatrix<REAL> &phiQ = data.phi;
  TPZFMatrix<REAL> &dphiQdaxes = data.dphix;//ELEMENTO DEFORMADO
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);

  TPZManVector<REAL,3> ax1(3),ax2(3), normal(3);
  for (int i=0; i<3; i++) {
    ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
    ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
  }
  Cross(ax1, ax2, normal);
  int phrq;
  phrq = data.fVecShapeIndex.NElements();
#ifdef LOG4CXX
	if(logger->isDebugEnabled() )
	{
		std::stringstream sout;
		sout<<std::endl;
		sout<<"el:"<<elId<<std::endl;
		sout<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  for(int iq=0; iq<phrq; iq++)
  {
    //ef(iq, 0) += 0.;
    int ivecind = data.fVecShapeIndex[iq].first;
    int ishapeind = data.fVecShapeIndex[iq].second;
    
    TPZManVector<REAL,3> ivecHDiv(3), ivecHCurl(3);
    for(int id=0; id<3; id++){
      ivecHDiv[id] = data.fNormalVec(id,ivecind);//JA EM XYZ
    }
    //ROTATE FOR HCURL
    Cross(normal, ivecHDiv, ivecHCurl);
		
    STATE ff = 0.;
    for (int i=0; i<3; i++) {
      ff += ivecHCurl[i]*force[i];
    }
//		if(elId==1){
//			ef(iq,0) += weight*ff*phiQ(ishapeind,0);
//		}
//		else{
//			ef(iq,0) -= weight*ff*phiQ(ishapeind,0);
//		}
    ef(iq,0) += weight*ff*phiQ(ishapeind,0);
    TPZManVector<REAL,3> curlI(3), gradPhiI(3);
    for (int i = 0; i<dphiQ.Rows(); i++) {
      gradPhiI[i] = dphiQ(i,ishapeind);
    }
    Cross(gradPhiI, ivecHCurl, curlI);
      STATE curlIStar;
    //eh necessario forcar curlE dot x = j*k0*sin(theta)*Ez
    curlIStar = -1. * imaginary * k0 * sin(fTheta) * phiQ(ishapeind,0) * ivecHCurl[2] ;
    
    for (int jq=0; jq<phrq; jq++)
    {
      TPZManVector<REAL,3> jvecHDiv(3), jvecHCurl(3);//JA EM XYZ
      int jvecind = data.fVecShapeIndex[jq].first;
      int jshapeind = data.fVecShapeIndex[jq].second;
      
      for(int id=0; id<3; id++){
        jvecHDiv[id] = data.fNormalVec(id,jvecind);
      }
      //ROTATE FOR HCURL
      Cross(normal, jvecHDiv,jvecHCurl);
      TPZManVector<REAL,3> curlJ(3),  gradPhiJ(3);
      
      for (int i = 0; i<dphiQ.Rows(); i++) {
        gradPhiJ[i] = dphiQ(i,jshapeind);
      }
      Cross(gradPhiJ, jvecHCurl, curlJ);
      //eh necessario forcar curlE dot x = j*k0*sin(theta)*Ez
        STATE curlJStar;

      curlJStar = imaginary * k0 * sin(fTheta) * phiQ(jshapeind,0) * jvecHCurl[2] *fScale;//AQUIFRAN
      
      REAL phiIdotphiJ;
      phiIdotphiJ = ( phiQ(ishapeind,0) * ivecHCurl[0] * phiQ(jshapeind,0) * jvecHCurl[0])
      + ( phiQ(ishapeind,0) * ivecHCurl[1] * phiQ(jshapeind,0) * jvecHCurl[1])
      + ( phiQ(ishapeind,0) * ivecHCurl[2] * phiQ(jshapeind,0) * jvecHCurl[2]);
      
        std::cout << " x " << data.x << std::endl;
        data.fNormalVec.Print("NormalVec");
        std::cout << " i " << iq << " j " << jq << " phiIdotphiJ " << phiIdotphiJ << std::endl;
        std::cout << " ivec " << ivecHCurl << " jvec " << jvecHCurl << std::endl;
        std::cout << " curli " << curlI << " curlj " << curlJ << std::endl;
        STATE stiff = 0.;
      stiff += (1./muR) * ( curlJStar * curlIStar + (curlJ[0] * curlI[0] + curlJ[1] * curlI[1] + curlJ[2] * curlI[2]) );
      stiff += -1. * k0 * k0 /(fScale * fScale) * epsilonR * phiIdotphiJ;
      
      ek(iq,jq) += weight * stiff;
    }
  }
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran2::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
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
        const STATE rhs = phiQ(i,0) * BIG * (1. + imaginary ) * v2;
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
          const STATE stiff = phiQ(i,0) * phiQ(j,0) * BIG * (1. + imaginary );
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
        STATE rhs = phiQ(i,0);
        rhs *= v2/fScale;
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
          STATE stiff = phiQ(i,0) *  phiQ(j,0);
          stiff *= v1/fScale;
          ek(i,j) += stiff*weight;
        }
      }
      break;
  }
}

void TPZMatValidacaoHCurlFran2::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatValidacaoHCurlFran2::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

int TPZMatValidacaoHCurlFran2::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2+elPMaxOrder*2;
}


int TPZMatValidacaoHCurlFran2::VariableIndex(const std::string &name)
{
  if(name == "absE") {
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
int TPZMatValidacaoHCurlFran2::NSolutionVariables(int var)
{
  int nVar = 0;
  switch (var) {
    case 2://absE
      nVar = 3;
      break;
    case 3://realE
      nVar = 3;
      break;
    default:
      nVar = TPZMaterial::NSolutionVariables(var);
      break;
  }
  return nVar;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatValidacaoHCurlFran2::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
  
  TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
  for (int i=0; i<3; i++) {
    ax1[i] = data.axes(0,i);
    ax2[i] = data.axes(1,i);
  }
  //ROTATE FOR HCURL
  Cross(ax1, ax2, normal);

  Solout.Resize(3);
  Cross(normal, data.sol[0], Solout);
  switch (var) {
    case 2://E
#ifdef STATE_COMPLEX
      Solout[0] = std::abs(Solout[0]);
      Solout[1] = std::abs(Solout[1]);
      Solout[2] = std::abs(Solout[2]);
#endif
      break;
    case 3://realE
#ifdef STATE_COMPLEX
      Solout[0] = std::real(Solout[0]);
      Solout[1] = std::real(Solout[1]);
      Solout[2] = std::real(Solout[2]);
#endif
      break;
    default:
      break;
  }
}


STATE urDefault( TPZVec<REAL> &x )
{
  return 1.0;
}

STATE erDefault( TPZVec<REAL> &x )
{
  return 1.0;
}









