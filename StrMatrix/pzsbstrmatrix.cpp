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
TPZMatrix<TVar> * TPZSBandStructMatrix<TVar,TPar>::CreateAssemble(TPZBaseMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<TVar> *mat = Create();
	rhs.Redim(mat->Rows(),1);
	Assemble(*mat,rhs,guiInterface);
    return mat;
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSBandStructMatrix<TVar,TPar>::Create(){
    if (fEquationFilter.IsActive()) {
        DebugStop();
    }
	int64_t neq = fEquationFilter.NActiveEquations();
	
	int64_t band = fMesh->BandWidth();
	return new TPZSBMatrix<TVar>(neq,band);
}

template<class TVar, class TPar>
TPZSBandStructMatrix<TVar,TPar>::TPZSBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

template<class TVar, class TPar>
TPZSBandStructMatrix<TVar,TPar>::TPZSBandStructMatrix(TPZAutoPointer<TPZCompMesh> mesh) : TPZStructMatrix(mesh)
{
}

template<class TVar, class TPar>
int TPZSBandStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSBandStructMatrix") ^
        TPZStructMatrix::ClassId() << 1 ^
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
