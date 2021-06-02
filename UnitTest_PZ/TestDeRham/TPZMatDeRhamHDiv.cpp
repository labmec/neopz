#include "TPZMatDeRhamHDiv.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"

void
TPZMatDeRhamHDiv::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef)
{
  
  const TPZFMatrix<REAL> &divPhi = data.divphi;
  const int nFuncs = divPhi.Rows();
  
  for (int i = 0; i < nFuncs; i++) {
    for (int j = 0; j < nFuncs; j++) {
      ek(i, j) += divPhi(i,0) * divPhi(j,0) * weight;
    }
  }
}