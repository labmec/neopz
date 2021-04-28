/**
 * @file
 * @brief Contains the implementation of the TPZSBandStructMatrix methods. 
 */

#include "pzsbstrmatrix.h"

#include "pzsbndmat.h"
#include "pzcmesh.h"
#include "TPZGuiInterface.h"

template<class TVar, class TPar>
TPZStructMatrix * TPZSBandStructMatrix<TVar,TPar>::Clone(){
    return new TPZSBandStructMatrix(*this);
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSBandStructMatrix<TVar,TPar>::Create(){
    if (this->fEquationFilter.IsActive()) {
        DebugStop();
    }
	int64_t neq = this->fEquationFilter.NActiveEquations();
	
	int64_t band = this->fMesh->BandWidth();
	return new TPZSBMatrix<TVar>(neq,band);
}

template<class TVar, class TPar>
int TPZSBandStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSBandStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSBandStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSBandStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSBandStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSBandStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSBandStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZSBandStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZSBandStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZSBandStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
