#include "TPZMatL2Product.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"

void TPZMatL2Product::Contribute(const TPZMaterialDataT<STATE> &data,
                                REAL weight, TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef){

  if (data.fShapeType == TPZMaterialData::EVecShape){
    ContributeVecShape(data, weight, ek, ef);
  }
  else if (data.fShapeType == TPZMaterialData::EScalarShape){
    ContributeScalarShape(data, weight, ek, ef);
  }
}

void TPZMatL2Product::ContributeVecShape(const TPZMaterialDataT<STATE> &data,
                                        REAL weight, TPZFMatrix<STATE> &ek,
                                        TPZFMatrix<STATE> &ef)
{
  int64_t nShapeHDiv = data.fVecShapeIndex.NElements(); // number of Hdiv shape functions
  TPZFNMatrix<150, REAL> phiHDiv(fDim, nShapeHDiv, 0.0);

  for (int i = 0; i < nShapeHDiv; i++){
    int ivecind = data.fVecShapeIndex[i].first;
    for (int j = 0; j < fDim; j++){
      phiHDiv(j, i) = data.fDeformedDirections(j, ivecind);
    }
  }

#ifdef PZ_USING_LAPACK
  ek.AddContribution(0, 0, phiHDiv, true, phiHDiv, false, weight);
#else
  for (int i = 0; i < nShapeHDiv; i++){
    for (int j = 0; j < nShapeHDiv; j++){
      STATE phiiphij = 0.0;
      for (int d = 0; d < fDim; d++){
        phiiphij += phiHDiv(d, i) * phiHDiv(d, j);
      }
      ek(i, j) += phiiphij * weight;
    }
  }
#endif
}

void TPZMatL2Product::ContributeScalarShape(const TPZMaterialDataT<STATE> &data,
                                           REAL weight, TPZFMatrix<STATE> &ek,
                                           TPZFMatrix<STATE> &ef){

  const TPZFMatrix<REAL> &phi = data.phi;
  const int nFuncs = phi.Rows();
  for (int i = 0; i < nFuncs; i++){
    for (int j = 0; j < nFuncs; j++){
      ek(i, j) += phi(i, 0) * phi(j, 0) * weight;
    }
  }
}