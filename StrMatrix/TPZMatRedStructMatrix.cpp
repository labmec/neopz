/**
 * @file
 * @brief Contains the implementation of TPZMatRedStructMatrix methods.
 * @author Giovane Avancini
 * @date 20/09/2023
 */

#include "TPZMatRedStructMatrix.h"
#include "pzcmesh.h"
#include "TPZMatRedMatrix.h"
#include "TPZRenumbering.h"
#include "TPZGuiInterface.h"
#include "TPZTimer.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template <class TVar, class TPar>
TPZStructMatrix *TPZMatRedStructMatrix<TVar, TPar>::Clone()
{
    return new TPZMatRedStructMatrix(*this);
}

template <class TVar, class TPar>
TPZMatrix<TVar> *TPZMatRedStructMatrix<TVar, TPar>::Create()
{
    //Needs implementation
    return nullptr;
}

template <class TVar, class TPar>
int TPZMatRedStructMatrix<TVar, TPar>::ClassId() const
{
    return Hash("TPZMatRedStructMatrix") ^
           TPZStructMatrixT<TVar>::ClassId() << 1 ^
           TPar::ClassId() << 2;
}

template <class TVar, class TPar>
void TPZMatRedStructMatrix<TVar, TPar>::Read(TPZStream &buf, void *context)
{
    TPZStructMatrix::Read(buf, context);
    TPar::Read(buf, context);
}

template <class TVar, class TPar>
void TPZMatRedStructMatrix<TVar, TPar>::Write(TPZStream &buf, int withclassid) const
{
    TPZStructMatrix::Write(buf, withclassid);
    TPar::Write(buf, withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZMatRedStructMatrix<STATE, TPZStructMatrixOR<STATE>>;
template class TPZMatRedStructMatrix<STATE, TPZStructMatrixOT<STATE>>;
template class TPZMatRedStructMatrix<STATE, TPZStructMatrixTBBFlow<STATE>>;

#ifndef USING_EIGEN
template class TPZMatRedStructMatrix<CSTATE, TPZStructMatrixOR<CSTATE>>;
template class TPZMatRedStructMatrix<CSTATE, TPZStructMatrixOT<CSTATE>>;
template class TPZMatRedStructMatrix<CSTATE, TPZStructMatrixTBBFlow<CSTATE>>;
#endif
