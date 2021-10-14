#include "TPZMatDeRhamH1HCurl.h"
#include "TPZCompElHCurl.h"
#include "pzaxestools.h"

TPZMatDeRhamH1HCurl* TPZMatDeRhamH1HCurl::NewMaterial() const{
  return new TPZMatDeRhamH1HCurl(*this);
}

void TPZMatDeRhamH1HCurl::Contribute(
  const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
  const TPZFMatrix<REAL> &gradPhiH1axes = datavec[fH1MeshIndex].dphix;
  TPZFNMatrix<3,REAL> gradPhiH1(3, gradPhiH1axes.Rows(), 0.);
  TPZAxesTools<REAL>::Axes2XYZ(gradPhiH1axes, gradPhiH1, datavec[fH1MeshIndex].axes);
  
  const auto &phiHCurl  = datavec[fHCurlMeshIndex].phi;

  const int nHCurl  = phiHCurl.Rows();
  const int nH1  = gradPhiH1.Cols();
  //position of first h1 func
  const int firstH1 = fH1MeshIndex * nHCurl;
  //position of first hcurl func
  const int firstHCurl = fHCurlMeshIndex * nH1;
    
  for(int iHCurl = 0; iHCurl < nHCurl; iHCurl++){
    for(int jHCurl = 0; jHCurl < nHCurl; jHCurl++){
      STATE phiIphiJ = 0;
      for(auto x = 0; x < fDim; x++){
        phiIphiJ += phiHCurl(iHCurl,x) * phiHCurl(jHCurl,x);
      }
      ek(firstHCurl+iHCurl,firstHCurl+jHCurl) += phiIphiJ * weight;
    }
    for(int jH1 = 0; jH1 < nH1; jH1++){
      STATE phiIgradPhiJ = 0;
      for(auto x = 0; x < fDim; x++){
        phiIgradPhiJ += phiHCurl(iHCurl,x) * gradPhiH1(x,jH1);
      }
      ek(firstHCurl+iHCurl,firstH1+jH1) += phiIgradPhiJ * weight;
    }
  }

  for (int iH1 = 0; iH1 < nH1; iH1++) {
    for (int jHCurl = 0; jHCurl < nHCurl; jHCurl++) {
      STATE gradPhiIphiJ = 0;
      for(auto x = 0; x < fDim; x++){
        gradPhiIphiJ += gradPhiH1(x,iH1) * phiHCurl(jHCurl,x);
      }
      ek(firstH1 + iH1, firstHCurl + jHCurl) += gradPhiIphiJ * weight;
    }
    for (int jH1 = 0; jH1 < nH1; jH1++) {
      STATE gradPhiIgradPhiJ = 0;
      for(auto x = 0; x < fDim; x++){
        gradPhiIgradPhiJ += gradPhiH1(x,iH1) * gradPhiH1(x,jH1);
      }
      ek(firstH1 + iH1, firstH1 + jH1) += gradPhiIgradPhiJ * weight;
    }
  }
}


int TPZMatDeRhamH1HCurl::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "SolutionLeft") == 0) return 0;
  if( strcmp(name.c_str(), "SolutionRight") == 0) return 1;
  DebugStop();
  return -1;
}

int TPZMatDeRhamH1HCurl::NSolutionVariables(int var) const
{
  switch (var) {
  case 0: // SolutionLeft
    return 1;
  case 1: // SolutionRight
    return fDim;
  default:
    DebugStop();
    break;
  }
  return 1;
}

void TPZMatDeRhamH1HCurl::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
              TPZVec<STATE> &solout)
{
  auto solH1 = datavec[fH1MeshIndex].sol[0];
  auto solHCurl = datavec[fHCurlMeshIndex].sol[0];

  switch (var) {
  case 0: {
    solout = solH1;
    break;
  }
  case 1: {
    solout = solHCurl;
    break;
  }
  default:
    DebugStop();
  }
}