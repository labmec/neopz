#include "TPZStrMatParInterface.h"
#include "pzbasematrix.h"
#include "TPZStructMatrix.h"
#include "TPZGuiInterface.h"

TPZBaseMatrix *TPZStrMatParInterface::CreateAssemble(
    TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    this->InitCreateAssemble();

    TPZStructMatrix *myself = dynamic_cast<TPZStructMatrix*>(this);
    TPZBaseMatrix *stiff = myself->Create();
    
    const int64_t cols = MAX(1, rhs.Cols());
    if(ComputeRhs()) rhs.Redim(myself->EquationFilter().NEqExpand(), cols);
    Assemble(*stiff, rhs, guiInterface);
    this->EndCreateAssemble(stiff);
#ifdef PZ_LOG2
    if (loggerel.isDebugEnabled()) {
        std::stringstream sout;
        stiff->Print("Stiffness matrix", sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel, sout.str())
    }
#endif
    return stiff;
}


int TPZStrMatParInterface::ClassId() const{
    return Hash("TPZStrMatParInterface");
}

void TPZStrMatParInterface::Read(TPZStream& buf, void* context) {
    buf.Read(&fNumThreads);
    buf.Read(fComputeRhs);
}

void TPZStrMatParInterface::Write(TPZStream& buf, int withclassid) const {
    buf.Write(&fNumThreads);
    buf.Write(fComputeRhs);
}