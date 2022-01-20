#include "TPZMatDeRhamHCurlHDiv.h"
#include "TPZCompElHCurl.h"
#include "pzaxestools.h"

TPZMatDeRhamHCurlHDiv* TPZMatDeRhamHCurlHDiv::NewMaterial() const{
  return new TPZMatDeRhamHCurlHDiv(*this);
}

void TPZMatDeRhamHCurlHDiv::Contribute(
  const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
  const TPZFMatrix<REAL> &deformed = datavec[fHDivMeshIndex].fDeformedDirections;
  
  const auto &curlphi  = datavec[fHCurlMeshIndex].curlphi;

  const int nHCurl  = curlphi.Cols();
  const int nHDiv  = deformed.Cols();
  //position of first h1 func
  const int firstHDiv = fHDivMeshIndex * nHCurl;
  //position of first hcurl func
  const int firstHCurl = fHCurlMeshIndex * nHDiv;
    
  for(int iHCurl = 0; iHCurl < nHCurl; iHCurl++){
    for(int jHCurl = 0; jHCurl < nHCurl; jHCurl++){
      STATE phiIphiJ = 0;
      for(auto x = 0; x < fDim; x++){
        phiIphiJ += curlphi(x,iHCurl) * curlphi(x,jHCurl);
      }
      ek(firstHCurl+iHCurl,firstHCurl+jHCurl) += phiIphiJ * weight;
    }
    for(int jHDiv = 0; jHDiv < nHDiv; jHDiv++){
      STATE phiIgradPhiJ = 0;
      for(auto x = 0; x < fDim; x++){
        phiIgradPhiJ += curlphi(x,iHCurl) * deformed(x,jHDiv);
      }
      ek(firstHCurl+iHCurl,firstHDiv+jHDiv) += phiIgradPhiJ * weight;
      ek(firstHDiv+jHDiv,firstHCurl+iHCurl) += phiIgradPhiJ * weight;
    }
  }

  for (int iHDiv = 0; iHDiv < nHDiv; iHDiv++) {
    for (int jHDiv = 0; jHDiv < nHDiv; jHDiv++) {
      STATE gradPhiIgradPhiJ = 0;
      for(auto x = 0; x < 3; x++){
        gradPhiIgradPhiJ += deformed(x,iHDiv) * deformed(x,jHDiv);
      }
      ek(firstHDiv + iHDiv, firstHDiv + jHDiv) += gradPhiIgradPhiJ * weight;
    }
  }
}


int TPZMatDeRhamHCurlHDiv::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "SolutionLeft") == 0) return 0;
  if( strcmp(name.c_str(), "SolutionRight") == 0) return 1;
  DebugStop();
  return -1;
}

int TPZMatDeRhamHCurlHDiv::NSolutionVariables(int var) const
{
  switch (var) {
  case 0: // SolutionLeft
    return 3;
  case 1: // SolutionRight
    return 3;
  default:
    DebugStop();
    break;
  }
  return 1;
}

void TPZMatDeRhamHCurlHDiv::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
              TPZVec<STATE> &solout)
{
  auto solH1 = datavec[fHDivMeshIndex].sol[0];
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