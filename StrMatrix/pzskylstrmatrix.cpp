/**
 * @file
 * @brief Contains the implementation of the TPZSkylineStructMatrix methods. 
 */

#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzcmesh.h"

TPZStructMatrix * TPZSkylineStructMatrix::Clone(){
    return new TPZSkylineStructMatrix(*this);
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

int TPZSkylineStructMatrix::ClassId() const{
    return Hash("TPZSkylineStructMatrix") ^ TPZStructMatrix::ClassId() << 1;
}

void TPZSkylineStructMatrix::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPZStructMatrixOR::Read(buf,context);
}

void TPZSkylineStructMatrix::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPZStructMatrixOR::Write(buf,withclassid);
}