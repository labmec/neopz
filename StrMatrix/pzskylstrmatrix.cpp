/**
 * @file
 * @brief Contains the implementation of the TPZSkylineStructMatrix methods. 
 */

#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzvec.h"

TPZStructMatrix * TPZSkylineStructMatrix::Clone(){
    return new TPZSkylineStructMatrix(*this);
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(const TPZSkylineStructMatrix &cp)
:TPZStructMatrix(cp) {
	//nothing here
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh)
{
}

TPZMatrix<STATE> * TPZSkylineStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	int neq = fEquationFilter.NEq();
	TPZMatrix<STATE> *stiff = Create();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs, guiInterface);
    return stiff;
}

TPZMatrix<STATE> * TPZSkylineStructMatrix::Create(){
    int neq = fMesh->NEquations();
    TPZVec<int> skyline;
	if (fOnlyInternal) {
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (fMesh);
		if (submesh) {
			submesh->SkylineInternal(skyline);
		}
		else {
			fMesh->Skyline(skyline);
		}
	}
	else {
		fMesh->Skyline(skyline);
	}
    fEquationFilter.FilterSkyline(skyline);
    neq = fEquationFilter.NEq();
	
    return this->ReallyCreate(neq,skyline);//new TPZSkylMatrix<STATE>(neq,skyline);
}


TPZMatrix<STATE> * TPZSkylineStructMatrix::ReallyCreate(int neq, const TPZVec<int> &skyline){
    return new TPZSkylMatrix<STATE>(neq,skyline);
}




TPZSkylineStructMatrix::~TPZSkylineStructMatrix(){
}




