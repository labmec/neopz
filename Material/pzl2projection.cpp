//$Id: pzl2projection.cpp,v 1.6 2007-07-04 19:04:28 tiago Exp $ 

#include "pzl2projection.h"

TPZL2Projection::TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol)
  :TPZDiscontinuousGalerkin(id){
  this->fDim = dim;
  this->fNStateVars = nstate;
  this->fSol = sol;
  this->SetIsReferred(false);
}

TPZL2Projection::TPZL2Projection(const TPZL2Projection &cp):TPZDiscontinuousGalerkin(cp){
  this->fDim = cp.fDim;
  this->fNStateVars = cp.fNStateVars;
  this->fSol = cp.fSol;
  this->SetIsReferred(cp.fIsReferred);
}

TPZL2Projection::~TPZL2Projection()
{
}

void TPZL2Projection::SetIsReferred(bool val){
  this->fIsReferred = val;
}

TPZAutoPointer<TPZMaterial> TPZL2Projection::NewMaterial(){
  return new TPZL2Projection(*this);
}

void TPZL2Projection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){

  if (this->HasForcingFunction()){
    this->fForcingFunction(data.x, this->fSol);
  }

  const int nvars = this->fNStateVars;
  if (this->fIsReferred){
    this->fSol.Resize(nvars);
    if (data.sol.NElements() < 2*nvars){//it means the referred element does not exist or it is an interface element. In that case, I ASSUME the referred solution is ZERO
      this->fSol.Fill(0.);
    }//if
    else{
      for(int i = 0; i < nvars; i++){
        this->fSol[i] = data.sol[i+nvars];
      }//for
    }//else
  }//if

  const int nshape = data.phi.Rows();
  for(int i = 0; i < nshape; i++){
    for(int j = 0; j < nshape; j++){
      for(int ivi = 0; ivi < nvars; ivi++){
        for(int ivj = 0; ivj < nvars; ivj++){
          const int posI = nvars*i+ivi;
          const int posJ = nvars*j+ivj;
          ek(posI, posJ) += weight*data.phi(i,0)*data.phi(j,0);
        }//ivj
      }//ivi
    }//for j
    for(int ivi = 0; ivi < nvars; ivi++){
      const int posI = nvars*i+ivi;
      ef(posI,0) += weight*data.phi(i,0)*this->fSol[ivi];
    }//ivi
  }//for i
}

