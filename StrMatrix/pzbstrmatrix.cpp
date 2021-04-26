/**
 * @file
 * @brief Contains the implementation of the TPZBandStructMatrix methods. 
 */

#include "pzbstrmatrix.h"
#include "pzbndmat.h"
#include "pzcmesh.h"
#include "TPZGuiInterface.h"

template<class TVar, class TPar>
TPZStructMatrix * TPZBandStructMatrix<TVar,TPar>::Clone(){
    return new TPZBandStructMatrix(*this);
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZBandStructMatrix<TVar,TPar>::Create(){
    if (this->fEquationFilter.IsActive()) {
        DebugStop();
    }
    int64_t neq = this->fEquationFilter.NActiveEquations();
    int64_t band = this->fMesh->BandWidth();
    return new TPZFBMatrix<TVar>(neq,band);
}

template<class TVar, class TPar>
int TPZBandStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZBandStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZBandStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZBandStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZBandStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZBandStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZBandStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZBandStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZBandStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZBandStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
