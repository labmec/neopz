//$Id: pzl2projection.cpp,v 1.4 2007-05-21 19:48:15 tiago Exp $ 

#include "pzl2projection.h"

TPZL2Projection::TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol)
  :TPZMaterial(id){
  this->fDim = dim;
  this->fNStateVars = nstate;
  this->fSol = sol;
  this->SetIsReferred(false);
}


TPZL2Projection::~TPZL2Projection()
{
}

void TPZL2Projection::SetIsReferred(bool val){
  this->fIsReferred = val;
}

void TPZL2Projection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){

  if (this->HasForcingFunction()){
    this->fForcingFunction(data.x, this->fSol);
  }

  const int nvars = this->fNStateVars;
  if (this->fIsReferred){
    this->fSol.Resize(nvars);
    for(int i = 0; i < nvars; i++){
      this->fSol[i] = data.sol[i+nvars];
    }//for
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

