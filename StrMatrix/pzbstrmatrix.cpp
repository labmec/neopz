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

TPZMatrix * TPZBandStructMatrix::CreateAssemble(TPZFMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix *stiff = Create();
	int neq = stiff->Rows();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs,guiInterface);
	return stiff;
}

TPZMatrix * TPZBandStructMatrix::Create(){
    int neq = fMesh->NEquations();
    if(HasRange()) neq = fMaxEq-fMinEq;
    int band = fMesh->BandWidth();
    return new TPZFBMatrix(neq,band);
}

TPZBandStructMatrix::TPZBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}
