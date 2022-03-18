#include "TPZMatDeRhamHDivL2.h"
#include "TPZCompElHCurl.h"
#include "pzaxestools.h"

TPZMatDeRhamHDivL2* TPZMatDeRhamHDivL2::NewMaterial() const{
  return new TPZMatDeRhamHDivL2(*this);
}

void TPZMatDeRhamHDivL2::Contribute(
  const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
  if(datavec[0].fShapeType != TPZShapeData::EVecShape) DebugStop();
  if(datavec[1].fShapeType != TPZShapeData::EScalarShape) DebugStop();

  const TPZFMatrix<REAL> &phiHDiv = datavec[fHDivMeshIndex].divphi;
  
  const auto &phiL2  = datavec[fL2MeshIndex].fPhi;

  const int64_t nHDiv  = phiHDiv.Rows();
  const int64_t nL2  = phiL2.Rows();
    
  for(int64_t iHDiv = 0; iHDiv < nHDiv; iHDiv++){
    for(int jHDiv = 0; jHDiv < nHDiv; jHDiv++){
        STATE phiIphiJ = phiHDiv(iHDiv,0) * phiHDiv(jHDiv,0);    
        ek(iHDiv,jHDiv) += phiIphiJ * weight;
    }
    for(int jL2 = 0; jL2 < nL2; jL2++){
      STATE phiIphiJ =  phiL2(jL2,0) * phiHDiv(iHDiv,0);
      ek(iHDiv,nHDiv+jL2) += phiIphiJ * weight;
      ek(nHDiv+jL2,iHDiv) += phiIphiJ * weight;
    }
  }

  for (int64_t iL2 = 0; iL2 < nL2; iL2++) {
    for (int64_t jL2 = 0; jL2 < nL2; jL2++) {
      STATE PhiIPhiJ = phiL2(iL2,0) * phiL2(jL2,0);
      ek(nHDiv + iL2, nHDiv + jL2) += PhiIPhiJ * weight;
    }
  }
}


int TPZMatDeRhamHDivL2::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "SolutionLeft") == 0) return 0;
  if( strcmp(name.c_str(), "SolutionRight") == 0) return 1;
  DebugStop();
  return -1;
}

int TPZMatDeRhamHDivL2::NSolutionVariables(int var) const
{
  switch (var) {
  case 0: // SolutionLeft
    return 3;
  case 1: // SolutionRight
    return 1;
  default:
    DebugStop();
    break;
  }
  return 1;
}

void TPZMatDeRhamHDivL2::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
              TPZVec<STATE> &solout)
{
  auto solHDiv = datavec[fHDivMeshIndex].sol[0];
  auto solL2 = datavec[fL2MeshIndex].sol[0];

  switch (var) {
  case 0: {
    solout = solHDiv;
    break;
  }
  case 1: {
    solout = solL2;
    break;
  }
  default:
    DebugStop();
  }
}