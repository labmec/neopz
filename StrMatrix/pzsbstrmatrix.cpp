/**
 * @file
 * @brief Contains the implementation of the TPZSBandStructMatrix methods. 
 */

#include "pzsbstrmatrix.h"

#include "pzsbndmat.h"
#include "pzcmesh.h"
#include "TPZGuiInterface.h"

TPZStructMatrix * TPZSBandStructMatrix::Clone(){
    return new TPZSBandStructMatrix(*this);
}

TPZMatrix<STATE> * TPZSBandStructMatrix::CreateAssemble(TPZBaseMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<STATE> *mat = Create();
	rhs.Redim(mat->Rows(),1);
	Assemble(*mat,rhs,guiInterface);
    return mat;
}

TPZMatrix<STATE> * TPZSBandStructMatrix::Create(){
    if (fEquationFilter.IsActive()) {
        DebugStop();
    }
	int64_t neq = fEquationFilter.NActiveEquations();
	
	int64_t band = fMesh->BandWidth();
	return new TPZSBMatrix<STATE>(neq,band);
}

TPZSBandStructMatrix::TPZSBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZSBandStructMatrix::TPZSBandStructMatrix(TPZAutoPointer<TPZCompMesh> mesh) : TPZStructMatrix(mesh)
{
}

int TPZSBandStructMatrix::ClassId() const{
    return Hash("TPZSBandStructMatrix") ^ TPZStructMatrix::ClassId() << 1;
}

void TPZSBandStructMatrix::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPZStructMatrixOR::Read(buf,context);
}

void TPZSBandStructMatrix::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPZStructMatrixOR::Write(buf,withclassid);
}