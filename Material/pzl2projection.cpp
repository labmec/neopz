//$Id: pzl2projection.cpp,v 1.10 2008-10-23 10:37:52 fortiago Exp $ 

#include "pzl2projection.h"
#include "pzbndcond.h"

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

void TPZL2Projection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc){

  const int nvars = this->fNStateVars;
  TPZFMatrix &phi = data.phi;
  const int phr = phi.Rows();
  int in, jn, iv;
  
  switch (bc.Type()){
  
    /// Dirichlet condition
    case 0 : {      
      for(iv = 0; iv < nvars; iv++){
        for(in = 0 ; in < phr; in++) {
          ef(nvars*in+iv,0) += TPZMaterial::gBigNumber * bc.Val2()(iv,0) * phi(in,0) * weight;      
          for (jn = 0 ; jn < phr; jn++) {
            ek(nvars*in+iv,nvars*jn+iv) += TPZMaterial::gBigNumber * phi(in,0) * phi(jn,0) * weight;
          }///jn
        }///in
      }///iv
      break;
    }
         
    default:{
      std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
    }
  }///switch

}

int TPZL2Projection::VariableIndex(const std::string &name){
  if(!strcmp("Solution",name.c_str())) return ESolution;
  return TPZMaterial::VariableIndex(name);
}

int TPZL2Projection::NSolutionVariables(int var){
  const int nvars = this->NStateVariables();
  if(var == ESolution) return nvars;

  return 0;
}

void TPZL2Projection::Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                               TPZFMatrix &axes, int var, TPZVec<REAL> &Solout){
  if (var == ESolution){
    Solout = Sol;
    return;
  }

  Solout.Resize(0);
}


