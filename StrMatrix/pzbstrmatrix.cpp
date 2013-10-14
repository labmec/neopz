/**
 * @file
 * @brief Contains the implementation of the TPZBandStructMatrix methods. 
 */

#include "pzbstrmatrix.h"
#include "pzbndmat.h"
#include "pzcmesh.h"

TPZBandStructMatrix::~TPZBandStructMatrix(){}

TPZStructMatrix * TPZBandStructMatrix::Clone(){
    return new TPZBandStructMatrix(*this);
}

TPZMatrix<STATE> * TPZBandStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<STATE> *stiff = Create();
	long neq = stiff->Rows();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs,guiInterface);
	return stiff;
}

TPZMatrix<STATE> * TPZBandStructMatrix::Create(){
    if (fEquationFilter.IsActive()) {
        DebugStop();
    }
    long neq = fEquationFilter.NActiveEquations();
    long band = fMesh->BandWidth();
    return new TPZFBMatrix<STATE>(neq,band);
}

TPZBandStructMatrix::TPZBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}
