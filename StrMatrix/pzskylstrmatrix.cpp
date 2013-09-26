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

TPZMatrix<STATE> * TPZSkylineStructMatrix::Create(){
    TPZVec<long> skyline;
    fMesh->Skyline(skyline);
    fEquationFilter.FilterSkyline(skyline);
    long neq = fEquationFilter.NActiveEquations();
    return this->ReallyCreate(neq,skyline);//new TPZSkylMatrix<STATE>(neq,skyline);
}


TPZMatrix<STATE> * TPZSkylineStructMatrix::ReallyCreate(long neq, const TPZVec<long> &skyline){
    return new TPZSkylMatrix<STATE>(neq,skyline);
}




TPZSkylineStructMatrix::~TPZSkylineStructMatrix(){
}




