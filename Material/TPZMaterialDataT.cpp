#include "TPZMaterialDataT.h"

template<class TVar>
TPZMaterialDataT<TVar>::TPZMaterialDataT() :
    TPZMaterialData(), sol(1),dsol(1){
}


template<class TVar>
void TPZMaterialDataT<TVar>::ComputeFunctionDivergence()
{
    std::cout << __PRETTY_FUNCTION__ << " should not be called anymore\n";
    // Getting test and basis functions
    //TPZFMatrix<REAL> dphi_s       = dphi; // Derivative For H1  test functions
    TPZShapeData *shapedata = this;
    int n_phi_v = shapedata->fSDVecShapeIndex.NElements();
#ifdef PZDEBUG
    if(divphi.Rows() < n_phi_v) DebugStop();
#endif
    REAL det_jac = detjac;

    int i_vec = 0;
    int i_phi_s = 0;
    
    for (int iq = 0; iq < n_phi_v; iq++)
    {
        i_vec = shapedata->fSDVecShapeIndex[iq].first;
        i_phi_s = shapedata->fSDVecShapeIndex[iq].second;
        divphi(iq,0) = 0.;

        int n_dir = shapedata->fDPhi.Rows();
        divphi(iq,0) = 0.;
        for (int k = 0; k < n_dir; k++) {
            divphi(iq,0) +=  shapedata->fDPhi(k,i_phi_s)*shapedata->fMasterDirections(k,i_vec)/detjac;
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
