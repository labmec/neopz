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

TPZSkylineStructMatrix::TPZSkylineStructMatrix() : TPZStructMatrix()
{
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh)
{
}

TPZMatrix<STATE> * TPZSkylineStructMatrix::Create(){
    TPZVec<int64_t> skyline;
    fMesh->Skyline(skyline);
    fEquationFilter.FilterSkyline(skyline);
    int64_t neq = fEquationFilter.NActiveEquations();
//    std::cout << skyline << std::endl;
    return this->ReallyCreate(neq,skyline);//new TPZSkylMatrix<STATE>(neq,skyline);
}


TPZMatrix<STATE> * TPZSkylineStructMatrix::ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline){
    return new TPZSkylMatrix<STATE>(neq,skyline);
}




TPZSkylineStructMatrix::~TPZSkylineStructMatrix(){
}




