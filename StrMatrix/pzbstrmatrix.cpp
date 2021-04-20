/**
 * @file
 * @brief Contains the implementation of the TPZBandStructMatrix methods. 
 */

#include "pzbstrmatrix.h"
#include "pzbndmat.h"
#include "pzcmesh.h"
#include "TPZGuiInterface.h"

TPZStructMatrix * TPZBandStructMatrix::Clone(){
    return new TPZBandStructMatrix(*this);
}

TPZMatrix<STATE> * TPZBandStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<STATE> *stiff = Create();
	int64_t neq = stiff->Rows();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs,guiInterface);
	return stiff;
}

TPZMatrix<STATE> * TPZBandStructMatrix::Create(){
    if (fEquationFilter.IsActive()) {
        DebugStop();
    }
    int64_t neq = fEquationFilter.NActiveEquations();
    int64_t band = fMesh->BandWidth();
    return new TPZFBMatrix<STATE>(neq,band);
}

TPZBandStructMatrix::TPZBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZBandStructMatrix::TPZBandStructMatrix(TPZAutoPointer<TPZCompMesh>mesh) : TPZStructMatrix(mesh)
{
}

int TPZBandStructMatrix::ClassId() const{
    return Hash("TPZBandStructMatrix") ^
        TPZStructMatrix::ClassId() << 1;
}

void TPZBandStructMatrix::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPZStructMatrixOR::Read(buf,context);
}

void TPZBandStructMatrix::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPZStructMatrixOR::Write(buf,withclassid);
}