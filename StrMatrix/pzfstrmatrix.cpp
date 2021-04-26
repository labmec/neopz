/**
 * @file
 * @brief Contains the implementation of the TPZFStructMatrix methods. 
 */

#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.tpzfstructmatrix");
static TPZLogger loggerel("pz.strmatrix.element");
#endif


using namespace std;
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZFStructMatrix<TVar,TPar>::Create(){
	int64_t neq = this->fEquationFilter.NActiveEquations();
    
	return new TPZFMatrix<TVar>(neq,neq,0.);
}

template<class TVar, class TPar>
TPZStructMatrix * TPZFStructMatrix<TVar,TPar>::Clone(){
    return new TPZFStructMatrix(*this);
}


template<class TVar, class TPar>
int TPZFStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZFStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZFStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZFStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZFStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZFStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZFStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZFStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZFStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZFStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;