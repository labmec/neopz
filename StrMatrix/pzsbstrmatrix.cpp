/**
 * @file
 * @brief Contains the implementation of the TPZSBandStructMatrix methods. 
 */

#include "pzsbstrmatrix.h"

#include "pzsbndmat.h"
#include "pzcmesh.h"

TPZStructMatrix * TPZSBandStructMatrix::Clone(){
    return new TPZSBandStructMatrix(*this);
}

TPZMatrix * TPZSBandStructMatrix::CreateAssemble(TPZFMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix *mat = Create();
	rhs.Redim(mat->Rows(),1);
	Assemble(*mat,rhs,guiInterface);
    return mat;
}

TPZMatrix * TPZSBandStructMatrix::Create(){
	int neq = fMesh->NEquations();
	if(HasRange())
	{
		neq = fMaxEq-fMinEq;
	}
	else
	{
		fMinEq = 0;
		fMaxEq = neq;
	}
	
	int band = fMesh->BandWidth();
	return new TPZSBMatrix(neq,band);
}

TPZSBandStructMatrix::TPZSBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}
