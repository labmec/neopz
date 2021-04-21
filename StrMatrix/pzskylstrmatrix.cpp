/**
 * @file
 * @brief Contains the implementation of the TPZSkylineStructMatrix methods. 
 */

#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzcmesh.h"

template<class TVar, class TPar>
TPZStructMatrix * TPZSkylineStructMatrix<TVar,TPar>::Clone(){
    return new TPZSkylineStructMatrix(*this);
}

template<class TVar, class TPar>
TPZSkylineStructMatrix<TVar,TPar>::TPZSkylineStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

template<class TVar, class TPar>
TPZSkylineStructMatrix<TVar,TPar>::TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh)
{
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSkylineStructMatrix<TVar,TPar>::Create(){
    TPZVec<int64_t> skyline;
    fMesh->Skyline(skyline);
    fEquationFilter.FilterSkyline(skyline);
    int64_t neq = fEquationFilter.NActiveEquations();
//    std::cout << skyline << std::endl;
    return this->ReallyCreate(neq,skyline);//new TPZSkylMatrix<TVar>(neq,skyline);
}


template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSkylineStructMatrix<TVar,TPar>::ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline){
    return new TPZSkylMatrix<TVar>(neq,skyline);
}

template<class TVar, class TPar>
int TPZSkylineStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSkylineStructMatrix") ^
        TPZStructMatrix::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSkylineStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSkylineStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSkylineStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSkylineStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
