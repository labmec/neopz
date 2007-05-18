//$Id: pzl2projection.cpp,v 1.3 2007-05-18 20:05:06 tiago Exp $ 

#include "pzl2projection.h"

TPZL2Projection::TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol)
  :TPZMaterial(id){
  this->fDim = dim;
  this->fNStateVars = nstate;
  this->fSol = sol;
}


TPZL2Projection::~TPZL2Projection()
{
}

void TPZL2Projection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){

  if (this->HasForcingFunction()){
    this->fForcingFunction(data.x, this->fSol);
  }

  const int nshape = data.phi.Rows();
  const int nvars = this->fNStateVars;
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

