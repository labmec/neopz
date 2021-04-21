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
TPZMatrix<TVar> * TPZBandStructMatrix<TVar,TPar>::CreateAssemble(TPZFMatrix<TVar> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<TVar> *stiff = Create();
	int64_t neq = stiff->Rows();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs,guiInterface);
	return stiff;
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZBandStructMatrix<TVar,TPar>::Create(){
    if (fEquationFilter.IsActive()) {
        DebugStop();
    }
    int64_t neq = fEquationFilter.NActiveEquations();
    int64_t band = fMesh->BandWidth();
    return new TPZFBMatrix<TVar>(neq,band);
}
template<class TVar, class TPar>
TPZBandStructMatrix<TVar,TPar>::TPZBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

template<class TVar, class TPar>
TPZBandStructMatrix<TVar,TPar>::TPZBandStructMatrix(TPZAutoPointer<TPZCompMesh>mesh) : TPZStructMatrix(mesh)
{
}

template<class TVar, class TPar>
int TPZBandStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZBandStructMatrix") ^
        TPZStructMatrix::ClassId() << 1 ^
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
