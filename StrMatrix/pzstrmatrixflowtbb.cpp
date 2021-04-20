/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixTBBFlow methods.
 */

#include "pzstrmatrixflowtbb.h"

#include "TPZStructMatrixTBBFlowUtils.h"
#include "TPZMaterial.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger loggerel("pz.strmatrix.element");
#endif



#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


#include "run_stats_table.h"

static RunStatsTable stat_ass_graph_tbb("-ass_graph_tbb", "Run statistics table for the graph creation, coloring and tbb::flow::graph TPZStructMatrixTBBFlow.");
#ifndef USING_TBB
#define NO_TBB \
    PZError<<"The class TPZStructMatrixTBBFlow depends on the TBB library.\n";\
    PZError<<"Please reconfigure the NeoPZ library using:\n";\
    PZError<<"USING_TBB=ON"<<std::endl;\
    DebugStop();
#endif
TPZStructMatrixTBBFlow::TPZStructMatrixTBBFlow() : TPZStructMatrixBase(){
#ifndef USING_TBB
    NO_TBB
#endif
}
TPZStructMatrixTBBFlow::TPZStructMatrixTBBFlow(TPZCompMesh *mesh) : fMesh(mesh), fEquationFilter(mesh->NEquations()) {
#ifndef USING_TBB
    NO_TBB
#else
    fMesh = mesh;
    this->SetNumThreads(0);
    this->fFlowGraph = new TPZFlowGraph(this);
#endif
}

TPZStructMatrixTBBFlow::TPZStructMatrixTBBFlow(TPZAutoPointer<TPZCompMesh> cmesh) : fCompMesh(cmesh), fEquationFilter(cmesh->NEquations()) {
#ifndef USING_TBB
    NO_TBB
#else
    fMesh = cmesh.operator->();
    this->SetNumThreads(0);
    this->fFlowGraph = new TPZFlowGraph(this);
#endif
}

TPZStructMatrixTBBFlow::TPZStructMatrixTBBFlow(const TPZStructMatrixTBBFlow &copy) : fMesh(copy.fMesh), fEquationFilter(copy.fEquationFilter)
{

#ifndef USING_TBB
    NO_TBB
#else
        if (copy.fCompMesh) {
        fCompMesh = copy.fCompMesh;
    }
    fMaterialIds = copy.fMaterialIds;
    fNumThreads = copy.fNumThreads;    
    fFlowGraph = new TPZFlowGraph(*copy.fFlowGraph);
#endif
}



TPZStructMatrixTBBFlow::~TPZStructMatrixTBBFlow()
{
    if (fFlowGraph) {
        delete fFlowGraph;
    }
}

TPZBaseMatrix *TPZStructMatrixTBBFlow::Create() {
    std::cout << "TPZStructMatrixTBBFlow::Create should never be called\n";
    return 0;
}

TPZStructMatrixTBBFlow *TPZStructMatrixTBBFlow::Clone() {
    std::cout << "TPZStructMatrixTBBFlow::Clone should never be called\n";
    return 0;
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixTBBFlow::Assemble(TPZBaseMatrix & stiffness, TPZBaseMatrix & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense,rhs.Cols(),0.);
        this->MultiThread_Assemble(stiffness,rhsloc,guiInterface);
        fEquationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        this->MultiThread_Assemble(stiffness,rhs,guiInterface);
        
    }
    ass_stiff.stop();
}

void TPZStructMatrixTBBFlow::Assemble(TPZBaseMatrix & rhs_base,TPZAutoPointer<TPZGuiInterface> guiInterface){
    if(!dynamic_cast<TPZFMatrix<STATE>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" Incompatible types. Aborting...\n";
        DebugStop();
    }
    auto &rhs = dynamic_cast<TPZFMatrix<STATE> &>(rhs_base);
    ass_rhs.start();
    if(fEquationFilter.IsActive())
    {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
        int64_t neqexpand = fEquationFilter.NEqExpand();
        //TODONORM
        if(rhs.Rows() != neqexpand || Norm(rhs) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense,1,0.);
        this->MultiThread_Assemble(rhsloc,guiInterface);
        fEquationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        this->MultiThread_Assemble(rhs,guiInterface);
    }
    ass_rhs.stop();
}

TPZBaseMatrix * TPZStructMatrixTBBFlow::CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZBaseMatrix *stiff = Create();
    
    int64_t cols = MAX(1, rhs.Cols());
    rhs.Redim(fEquationFilter.NEqExpand(),cols);
    
    Assemble(*stiff,rhs,guiInterface);
    
#ifdef PZ_LOG2
    if(loggerel.isDebugEnabled())
    {
        std::stringstream sout;
        stiff->Print("Stiffness matrix",sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
#endif
    return stiff;
    
}

void TPZStructMatrixTBBFlow::MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    this->fFlowGraph->ExecuteGraph(&rhs, &mat);
#else
    std::cout << "To use the tbb flow graph assemble please compile the NeoPZ with USING_TBB." << std::endl;
    DebugStop();
#endif
}


void TPZStructMatrixTBBFlow::MultiThread_Assemble(TPZBaseMatrix & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    this->fFlowGraph->ExecuteGraph(&rhs);
#else
    std::cout << "To use the tbb flow graph assemble please compile the NeoPZ with USING_TBB." << std::endl;
    DebugStop();
#endif
}

int TPZStructMatrixTBBFlow::ClassId() const{
    return Hash("TPZStructMatrixTBBFlow") ^ TPZStructMatrixBase::ClassId() << 1;
}


void TPZStructMatrixTBBFlow::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf,context);
    fMesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
    fCompMesh = TPZAutoPointerDynamicCast<TPZCompMesh>(TPZPersistenceManager::GetAutoPointer(&buf));
    fEquationFilter.Read(buf, context);
    buf.Read(fMaterialIds);
    buf.Read(&fNumThreads);
    PZError<<__PRETTY_FUNCTION__<<" not implemented. Aborting..."<<std::endl;
    DebugStop();
	//fFlowGraph = new TPZFlowGraph(this);
}

void TPZStructMatrixTBBFlow::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf,withclassid);
    TPZPersistenceManager::WritePointer(fMesh, &buf);
    TPZPersistenceManager::WritePointer(fCompMesh.operator ->(), &buf);
    fEquationFilter.Write(buf, withclassid);
    PZError<<__PRETTY_FUNCTION__<<" not implemented. Aborting..."<<std::endl;
    DebugStop();
//    TPZPersistenceManager::WritePointer(fFlowGraph, &buf);
    buf.Write(fMaterialIds);
    buf.Write(&fNumThreads);
}

template class TPZRestoreClass<TPZStructMatrixTBBFlow>;
