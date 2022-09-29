#include "TPZL2ProjectionHDiv.h"

template<class TVar>
TPZL2ProjectionHDiv<TVar>::TPZL2ProjectionHDiv(int id, int dim, int nstate) : TPZL2Projection<TVar>(id,dim,nstate){

}

template<class TVar>
TPZL2ProjectionHDiv<TVar>::TPZL2ProjectionHDiv(int id, int dim, int nstate, const TPZVec<TVar> &sol) : TPZL2Projection<TVar>(id,dim,nstate,sol){

}

template<class TVar>
void TPZL2ProjectionHDiv<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                            REAL weight,
                                            TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
	const int nshape = data.phi.Rows();
    const int nvars = this->fNStateVars;
    TPZManVector<TVar,10> solLoc(this->fSol);
    if(this->HasForcingFunction()){
        this->fForcingFunction(data.x,solLoc);
    }
    const auto &phi = data.phi;

    int phrq = data.fVecShapeIndex.NElements();
    const STATE inv_perm = 1.;
	//Calculate the matrix contribution for flux. Matrix A
    for (int iq = 0; iq < phrq; iq++) {
        //ef(iq, 0) += 0.;
        int ivecind = data.fVecShapeIndex[iq].first;
        int ishapeind = data.fVecShapeIndex[iq].second;
        TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
        for (int id = 0; id < 3; id++) {
            ivec(id, 0) = data.fDeformedDirections(id, ivecind);
        }

        TPZFNMatrix<3, REAL> ivecZ(3, 1, 0.);
        TPZFNMatrix<3, REAL> jvecZ(3, 1, 0.);
        for (int jq = 0; jq < phrq; jq++) {
            TPZFNMatrix<3, REAL> jvec(3, 1, 0.);
            int jvecind = data.fVecShapeIndex[jq].first;
            int jshapeind = data.fVecShapeIndex[jq].second;

            for (int id = 0; id < 3; id++) {
                jvec(id, 0) = data.fDeformedDirections(id, jvecind);
            }

            //dot product between Kinv[u]v
            jvecZ.Zero();
            for (int id = 0; id < 3; id++) {
                jvecZ(id, 0) += inv_perm * jvec(id, 0);
            }
            REAL prod1 = ivec(0, 0) * jvecZ(0, 0) + ivec(1, 0) * jvecZ(1, 0) + ivec(2, 0) * jvecZ(2, 0);
            ek(iq, jq) += weight * phi(ishapeind, 0) * phi(jshapeind, 0) * prod1;
        }
    }

    for (int iq = 0; iq < phrq; iq++) {
        for (int idim = 0; idim < this->fDim; idim++) {
            ef(iq, 0) +=  solLoc[0] * data.fDeformedDirections(idim, iq) * weight;
        }
    }


}


template<class TVar>
int TPZL2ProjectionHDiv<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return TPZL2Projection<TVar>::ESolution;
	return TPZL2Projection<TVar>::VariableIndex(name);
}

template<class TVar>
int TPZL2ProjectionHDiv<TVar>::NSolutionVariables(int var) const{
	if(var == TPZL2Projection<TVar>::ESolution) return this->fDim;

    return TPZL2Projection<TVar>::NSolutionVariables(var);
}



template class TPZL2ProjectionHDiv<STATE>;