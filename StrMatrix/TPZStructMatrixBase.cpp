#include "TPZStructMatrixBase.h"

#include "pzcmesh.h"
#include "pzcheckmesh.h"
#include "pzerror.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "TPZThreadPool.h"

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixOR"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
static LoggerPtr loggerGlobStiff(Logger::getLogger("pz.strmatrix.globalstiffness"));
#endif

TPZStructMatrixBase::TPZStructMatrixBase() : fMesh(NULL), fEquationFilter(0) {
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

TPZStructMatrixBase::TPZStructMatrixBase(TPZCompMesh *mesh)
    : fEquationFilter(0) {
    SetMesh(mesh);
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

TPZStructMatrixBase::TPZStructMatrixBase(TPZAutoPointer<TPZCompMesh> mesh)
    : fEquationFilter(0) {
    SetMesh(mesh);
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

TPZStructMatrixBase::TPZStructMatrixBase(const TPZStructMatrixBase &copy)
    : fMesh(copy.fMesh), fCompMesh(copy.fCompMesh),
      fEquationFilter(copy.fEquationFilter), fMaterialIds(copy.fMaterialIds),
      fNumThreads(copy.fNumThreads) {
}

void TPZStructMatrixBase::SetMesh(TPZCompMesh *mesh) {
    fMesh = mesh;
    fEquationFilter.SetNumEq(mesh ? mesh->NEquations() : 0);
#ifdef PZDEBUG
    if (mesh) {
        TPZCheckMesh checkmesh(fMesh, &std::cout);
        if (checkmesh.CheckConnectSeqNumberConsistency() != 0) {
            DebugStop();
        }
    }
#endif
}

void TPZStructMatrixBase::SetMesh(TPZAutoPointer<TPZCompMesh> mesh) {
    fCompMesh = mesh;
    SetMesh(mesh.operator->());
}

TPZMatrix<STATE> *TPZStructMatrixBase::CreateAssemble(
    TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    TPZMatrix<STATE> *stiff = Create();
    
    int64_t cols = MAX(1, rhs.Cols());
    rhs.Redim(fEquationFilter.NEqExpand(), cols);
    Assemble(*stiff, rhs, guiInterface);

#ifdef LOG4CXX2
    if (loggerel->isDebugEnabled()) {
        std::stringstream sout;
        stiff->Print("Stiffness matrix", sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel, sout.str())
    }
#endif
    return stiff;
}

void TPZStructMatrixBase::FilterEquations(TPZVec<int64_t> &origindex, TPZVec<int64_t> &destindex) const
{
    //destindex = origindex;
    fEquationFilter.Filter(origindex, destindex);
    
}

void TPZStructMatrixBase::SetMaterialIds(const std::set<int> &materialids)
{
    fMaterialIds = materialids;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::set<int>::const_iterator it;
        std::stringstream sout;
        sout << "setting input material ids ";
        for(it=materialids.begin(); it!= materialids.end(); it++)
        {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    if(!fMesh)
    {
        LOGPZ_WARN(logger,"SetMaterialIds called without mesh")
        return;
    }
    int64_t iel;
    TPZAdmChunkVector<TPZCompEl*> &elvec = fMesh->ElementVec();
    int64_t nel = elvec.NElements();
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = elvec[iel];
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
        if(!subcmesh) continue;
        TPZAutoPointer<TPZAnalysis> anal = subcmesh->Analysis();
        if(!anal)
        {
            LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without analysis object")
            DebugStop();
        }
        TPZStructMatrixBase *str = anal->StructMatrix().operator->();
        if(!str)
        {
            LOGPZ_WARN(logger,"SetMaterialIds called for substructure without structural matrix")
            continue;
        }
        str->SetMaterialIds(materialids);
    }
}

int TPZStructMatrixBase::ClassId() const{
    return Hash("TPZStructMatrixBase");
}

void TPZStructMatrixBase::Read(TPZStream& buf, void* context) {
    fMesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
    fCompMesh = TPZAutoPointerDynamicCast<TPZCompMesh>(TPZPersistenceManager::GetAutoPointer(&buf));
    fEquationFilter.Read(buf,context);
    buf.Read(fMaterialIds);
    buf.Read(&fNumThreads);
}

void TPZStructMatrixBase::Write(TPZStream& buf, int withclassid) const {
    TPZPersistenceManager::WritePointer(fMesh, &buf);
    TPZPersistenceManager::WritePointer(fCompMesh.operator ->(), &buf);
    fEquationFilter.Write(buf,withclassid);
    buf.Write(fMaterialIds);
    buf.Write(&fNumThreads);
}
