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
  const TPZFMatrix<REAL> &phiH1 = datavec[fH1MeshIndex].phi;
  const TPZFMatrix<REAL> &gradPhiH1axes = datavec[fH1MeshIndex].dphix;
  TPZFNMatrix<3,REAL> gradPhiH1(3, phiH1.Rows(), 0.);
  TPZAxesTools<REAL>::Axes2XYZ(gradPhiH1axes, gradPhiH1, datavec[fH1MeshIndex].axes);
  
  TPZFNMatrix<30,REAL> phiHCurl;
  
  TPZHCurlAuxClass::ComputeShape(datavec[fHCurlMeshIndex].fVecShapeIndex,
                                 phiH1,
                                 datavec[fHCurlMeshIndex].fDeformedDirections,
                                 phiHCurl);

  const int nHCurl  = phiHCurl.Rows();
  const int nH1  = phiH1.Rows();
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