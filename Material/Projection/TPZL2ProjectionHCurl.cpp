#include "TPZL2ProjectionHCurl.h"

/*
  @fran: TODO: these materials are inconsistent.
  They must be rewritten taken into account that
  HCurl basis functions are always 3d and their curl
  is a scalar in 1d/2d
**/

template<class TVar>
TPZL2ProjectionHCurl<TVar>::TPZL2ProjectionHCurl(int id, int dim, int nstate) : TPZL2Projection<TVar>(id,dim,nstate){

}

template<class TVar>
TPZL2ProjectionHCurl<TVar>::TPZL2ProjectionHCurl(int id, int dim, int nstate, const TPZVec<TVar> &sol) : TPZL2Projection<TVar>(id,dim,nstate,sol){

}

template<class TVar>
void TPZL2ProjectionHCurl<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                            REAL weight,
                                            TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
	const int nshape = data.phi.Rows();
    const int nvars = this->fNStateVars;
    TPZManVector<TVar,10> solLoc(this->fSol);
    if(this->HasForcingFunction()){
        this->fForcingFunction(data.x,solLoc);
    }
    const auto &phi = data.phi;
    int phrq = nshape;
	
    for (int iq = 0; iq < phrq; iq++) {
        for (int jq = 0; jq < phrq; jq++) {
            for (int idim = 0; idim < this->fDim; idim++) {
                ek(iq, jq) += weight * ( phi(iq, idim) * phi(jq, idim) );
            }
        }
    }

    for (int iq = 0; iq < phrq; iq++) {
        for (int idim = 0; idim < this->fDim; idim++) {
            ef(iq, 0) +=  solLoc[0] * phi(iq, idim) * weight;
        }
    }
}


template<class TVar>
int TPZL2ProjectionHCurl<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return TPZL2Projection<TVar>::ESolution;
	return TPZL2Projection<TVar>::VariableIndex(name);
}

template<class TVar>
int TPZL2ProjectionHCurl<TVar>::NSolutionVariables(int var) const{
	if(var == TPZL2Projection<TVar>::ESolution) return this->fDim;

    return TPZL2Projection<TVar>::NSolutionVariables(var);
}



template class TPZL2ProjectionHCurl<STATE>;
