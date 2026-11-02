#include "TPZAlgebraicInterface.h"
#include "TPZElementMatrixT.h"
#include "pzgeoel.h"
#include "pzcmesh.h"

TPZAlgebraicInterface::TPZAlgebraicInterface() : TPZRegisterClassId(&TPZAlgebraicInterface::ClassId), TPZCompEl(), fConnectIndexes()
{
}

TPZAlgebraicInterface::TPZAlgebraicInterface(TPZCompMesh &mesh, TPZGeoEl *reference) : TPZRegisterClassId(&TPZAlgebraicInterface::ClassId), TPZCompEl(mesh, reference)
{
    fConnectIndexes.fill(-1);
}

TPZAlgebraicInterface::TPZAlgebraicInterface(TPZCompMesh &mesh, const TPZCompEl &copy) : TPZRegisterClassId(&TPZAlgebraicInterface::ClassId), TPZCompEl(mesh, copy)
{
    const TPZAlgebraicInterface *copyPtr = dynamic_cast<const TPZAlgebraicInterface *>(&copy);
    if (!copyPtr) {
        DebugStop();
    }
    fConnectIndexes = copyPtr->fConnectIndexes;
    fMultiplier = copyPtr->fMultiplier;
}

TPZAlgebraicInterface::TPZAlgebraicInterface(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int64_t,int64_t> &gl2lcElMap) :
    TPZRegisterClassId(&TPZAlgebraicInterface::ClassId), TPZCompEl(mesh, copy, gl2lcElMap), fConnectIndexes()
{
    const TPZAlgebraicInterface *copyPtr = dynamic_cast<const TPZAlgebraicInterface *>(&copy);
    if (!copyPtr) {
        DebugStop();
    }
    fConnectIndexes = copyPtr->fConnectIndexes;
    fMultiplier = copyPtr->fMultiplier;
}

TPZAlgebraicInterface::~TPZAlgebraicInterface()
{
}

TPZCompEl *TPZAlgebraicInterface::Clone(TPZCompMesh &mesh) const
{
    return new TPZAlgebraicInterface(mesh, *this);
}

TPZCompEl *TPZAlgebraicInterface::ClonePatchEl(TPZCompMesh &mesh,
                                               std::map<int64_t,int64_t> &gl2lcConMap,
                                               std::map<int64_t,int64_t> &gl2lcElMap) const
{
    TPZAlgebraicInterface *clone = new TPZAlgebraicInterface(mesh, *this, gl2lcElMap);
    for (auto &connect : clone->fConnectIndexes) {
        if (connect >= 0) {
            auto it = gl2lcConMap.find(connect);
            if (it != gl2lcConMap.end()) {
                connect = it->second;
            }
        }
    }
    return clone;
}

int TPZAlgebraicInterface::NConnects() const
{
    return 2;
}

int64_t TPZAlgebraicInterface::ConnectIndex(int i) const
{
    if (i < 0 || i >= 2) {
        DebugStop();
    }
    return fConnectIndexes[i];
}

int TPZAlgebraicInterface::Dimension() const
{
    
    return 0;
}

void TPZAlgebraicInterface::SetConnectIndex(int inode, int64_t index)
{
    if (inode < 0 || inode >= 2) {
        DebugStop();
    }
    fConnectIndexes[inode] = index;
}

void TPZAlgebraicInterface::CalcStiff(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef)
{
    CalcStiffInternal(ek, ef);
}

void TPZAlgebraicInterface::CalcStiff(TPZElementMatrixT<CSTATE> &ek, TPZElementMatrixT<CSTATE> &ef)
{
    CalcStiffInternal(ek, ef);
}

int TPZAlgebraicInterface::ClassId() const
{
    return Hash("TPZAlgebraicInterface") ^ TPZCompEl::ClassId() << 1;
}

template<class TVar>
void TPZAlgebraicInterface::CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef)
{
    ek.Reset(Mesh(),TPZElementMatrixT<TVar>::EK);
    ef.Reset(Mesh(),TPZElementMatrixT<TVar>::EF);
    TPZCompEl::InitializeElementMatrix(ek, ef);
    ek.fConnect.Resize(2);
    ef.fConnect.Resize(2);
    for (int i = 0; i < 2; i++) {
        if (fConnectIndexes[i] >= 0) {
            ek.fConnect[i] = fConnectIndexes[i];
            ef.fConnect[i] = fConnectIndexes[i];
        } else {
            DebugStop();
        }
    }
    int eksize = 0;
    int dofsize[2];
    for (int i = 0; i < 2; i++) {
        int64_t connectindex = fConnectIndexes[i];
        TPZConnect &c = Mesh()->ConnectVec()[connectindex];
        eksize += c.NDof();
        dofsize[i] = c.NDof();
    }
    if(dofsize[0] != dofsize[1]) {
        DebugStop();
    }
    ek.fMat.Redim(eksize, eksize);
    ef.fMat.Redim(eksize, 1);
    for(int i = 0; i < dofsize[0]; i++) {
        ek.fMat(i, i+dofsize[0]) = fMultiplier;
        ek.fMat(i+dofsize[0], i) = fMultiplier;
    }
}
