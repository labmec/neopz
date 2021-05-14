#include "TPZMaterialDataT.h"

template<class TVar>
TPZMaterialDataT<TVar>::TPZMaterialDataT() :
    TPZMaterialData(), sol(1),dsol(1){
}


template<class TVar>
void TPZMaterialDataT<TVar>::ComputeFunctionDivergence()
{
    
    // Getting test and basis functions
    TPZFMatrix<REAL> dphi_s       = dphi; // Derivative For H1  test functions
    
    int n_phi_v = fVecShapeIndex.NElements();
#ifdef PZDEBUG
    if(divphi.Rows() < n_phi_v) DebugStop();
#endif
    REAL det_jac = detjac;

    int i_vec = 0;
    int i_phi_s = 0;
    
    for (int iq = 0; iq < n_phi_v; iq++)
    {
        i_vec = fVecShapeIndex[iq].first;
        i_phi_s = fVecShapeIndex[iq].second;
        divphi(iq,0) = 0.;

        int n_dir = dphi_s.Rows();
        divphi(iq,0) = 0.;
        for (int k = 0; k < n_dir; k++) {
            divphi(iq,0) +=  dphi(k,i_phi_s)*fMasterDirections(k,i_vec)/detjac;
        }
    }
        

}

template<class TVar>
void TPZMaterialDataT<TVar>::SetSolSizes(const int nSol, const int uLen,
                                         const int duRow, const int duCol)
{
    this->sol.Resize(nSol);
    this->dsol.Resize(nSol);
    for( auto& u : this->sol){
        u.Resize(uLen);
    }
    for( auto& du : this->dsol){
        du.Resize(duRow,duCol);
    }
}

template class TPZMaterialDataT<STATE>;
template class TPZMaterialDataT<CSTATE>;