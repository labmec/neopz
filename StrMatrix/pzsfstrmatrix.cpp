/**
 * @file
 * @brief Contains the implementation of the TPZSFStructMatrix methods. 
 */

#include "pzsfstrmatrix.h"
#include "pzsfulmat.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.tpzfstructmatrix");
static TPZLogger loggerel("pz.strmatrix.element");
#endif


using namespace std;
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSFStructMatrix<TVar,TPar>::Create(){
	int64_t neq = this->fEquationFilter.NActiveEquations();
    
	return new TPZSFMatrix<TVar>(neq);
}

template<class TVar, class TPar>
TPZStructMatrix * TPZSFStructMatrix<TVar,TPar>::Clone(){
    return new TPZSFStructMatrix(*this);
}


template<class TVar, class TPar>
int TPZSFStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSFStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSFStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSFStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSFStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSFStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSFStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZSFStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZSFStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZSFStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;