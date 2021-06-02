#include "TPZMatDeRhamL2.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"

void TPZMatDeRhamL2::Contribute(const TPZMaterialDataT<STATE> &data,
                                REAL weight, TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef) {
  const TPZFMatrix<REAL> &phi = data.phi;
  const int nFuncs = phi.Rows();
  for (int i = 0; i < nFuncs; i++) {
    for (int j = 0; j < nFuncs; j++) {
      ek(i, j) += phi(i, 0) * phi(j, 0) * weight;
    }
  }
}