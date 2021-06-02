#include "TPZMatDeRhamHCurl.h"
#include "TPZMaterialDataT.h"
#include "TPZCompElHCurl.h"
#include "pzaxestools.h"

void
TPZMatDeRhamHCurl::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef)
{
  const TPZFMatrix<REAL> &curlPhi = data.curlphi;
  const int nFuncs = curlPhi.Cols();
  const int curlDim = 2*fDim - 3;
  
  for (int i = 0; i < nFuncs; i++) {
    for (int j = 0; j < nFuncs; j++) {
      for (int x = 0; x < curlDim; x++)
        ek(i, j) += curlPhi(x,i) * curlPhi(x,j) * weight;
    }
  }
}