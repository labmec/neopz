#include "TPZStructMatrixBase.h"

#include "pzcmesh.h"
#include "pzcheckmesh.h"
#include "pzerror.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "TPZThreadPool.h"
#include "pzshtmat.h"
#include "pzlog.h"
#ifdef PZ_LOG
static PZLogger logger("pz.strmatrix.TPZStructMatrixOR");
static PZLogger loggerel("pz.strmatrix.element");
static PZLogger loggerel2("pz.strmatrix.elementinterface");
static PZLogger loggerelmat("pz.strmatrix.elementmat");
static PZLogger loggerCheck("pz.strmatrix.checkconsistency");
static PZLogger loggerGlobStiff("pz.strmatrix.globalstiffness");
#endif

//this is not in header file to avoid including cmesh.h there
TPZStructMatrixBase::~TPZStructMatrixBase(){}

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

void TPZStructMatrixBase::FilterEquations(TPZVec<int64_t> &origindex, TPZVec<int64_t> &destindex) const
{
    //destindex = origindex;
    fEquationFilter.Filter(origindex, destindex);
    
}

void TPZStructMatrixBase::SetMaterialIds(const std::set<int> &materialids)
{
    fMaterialIds = materialids;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
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
        TPZAutoPointer<TPZAnalysis> analysis = subcmesh->Analysis();
        if(!analysis)
        {
            LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without analysis object")
            DebugStop();
        }
        TPZStructMatrixBase *str = analysis->StructMatrix().operator->();
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

/// compute a color for each element
// @return the number of colors for parallel assembly
// the color =-1 when the element should not be computed
int TPZStructMatrixBase::ComputeElementColors(TPZVec<int> &elementcolors)
{
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    if(!fMaterialIds.size())
    {
        fMesh->ComputeElGraph(elgraph, elgraphindex);
    }
    else
    {
        fMesh->ComputeElGraph(elgraph, elgraphindex, fMaterialIds);
    }
    int64_t nconnects = fMesh->NConnects();
    const int color_stored = 10;
    TPZGenMatrix<int64_t> domain(nconnects,color_stored);
    int64_t nel = elgraphindex.size()-1;
    elementcolors.resize(nel);
    elementcolors.Fill(-1);
    bool element_not_colored = true;
    int highest_color = 0;
    int min_color = 0;
    int max_color = min_color+color_stored;
    while(element_not_colored)
    {
        element_not_colored = false;
        for (int64_t el = 0; el<nel; el++) {
            int icolor;
            bool color_found = true;
            for(int icolor = 0; icolor<color_stored; icolor++)
            {
                color_found = true;
                int64_t firstnode = elgraphindex[el];
                int64_t lastnode = elgraphindex[el+1];
                for (int inod = firstnode; inod < lastnode; inod++) {
                    int64_t ic = elgraph[inod];
                    if(domain(ic,icolor) != 0)
                    {
                        color_found = false;
                        break;
                    }
                    domain(ic,icolor) = 1;
                }
                if(color_found == true)
                {
                    elementcolors[el] = min_color+icolor;
                    if(min_color+icolor > highest_color) highest_color = min_color+icolor;
                    break;
                }
            }
            // the element didn't fit any color, the loop will have to be executed again
            if(!color_found) element_not_colored = true;
        }
        min_color += color_stored;
        max_color += color_stored;
        domain.Fill(0);
    }
    return highest_color;
}

