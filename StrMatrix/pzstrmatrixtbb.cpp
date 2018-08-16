
#include "TPZPersistenceManager.h"

/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixTBB methods.
 */
#ifdef USING_TBB

#include "pzstrmatrixtbb.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "pzsfulmat.h"

#include "pzgnode.h"
#include "TPZTimer.h"
#include "TPZThreadTools.h"


#include "pzcheckconsistency.h"
#include "TPZMaterial.h"

#include "pzlog.h"

#include "pz_pthread.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixTBB"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


#include "run_stats_table.h"


static RunStatsTable stat_ass_graph_tbb("-ass_graph_tbb", "Run statistics table for the graph creation, coloring and tbb::flow::graph TPZStructMatrixTBB.");

TPZStructMatrixTBB::TPZStructMatrixTBB() : TPZStructMatrixBase() {
	this->SetNumThreads(0);
#ifdef USING_TBB
	this->fFlowGraph = NULL;
#endif
}

TPZStructMatrixTBB::TPZStructMatrixTBB(TPZCompMesh *mesh, bool onlyrhs) : TPZStructMatrixBase(mesh) {
    this->SetNumThreads(0);
	if (mesh)
	   fEquationFilter = mesh->NEquations();
#ifdef USING_TBB
    this->fFlowGraph = new TPZFlowGraph(this,onlyrhs);
#endif
}

TPZStructMatrixTBB::TPZStructMatrixTBB(TPZAutoPointer<TPZCompMesh> cmesh, bool onlyrhs) : TPZStructMatrixBase(cmesh) {
    this->SetNumThreads(0);
#ifdef USING_TBB
    this->fFlowGraph = new TPZFlowGraph(this, onlyrhs);
#endif
}

TPZStructMatrixTBB::TPZStructMatrixTBB(const TPZStructMatrixTBB &copy) : TPZStructMatrixBase(copy)
{
#ifdef USING_TBB
    fFlowGraph = new TPZFlowGraph(*copy.fFlowGraph);
#endif
}



TPZStructMatrixTBB::~TPZStructMatrixTBB()
{
#ifdef USING_TBB
    if (fFlowGraph) {
        delete fFlowGraph;
    }
#endif
}

TPZMatrix<STATE> *TPZStructMatrixTBB::Create() {
    std::cout << "TPZStructMatrixTBB::Create should never be called\n";
    return 0;
}

TPZStructMatrixTBB *TPZStructMatrixTBB::Clone() {
    return new TPZStructMatrixTBB(*this);
    
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixTBB::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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

void TPZStructMatrixTBB::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_rhs.start();
    if(fEquationFilter.IsActive())
    {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
        int64_t neqexpand = fEquationFilter.NEqExpand();
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

TPZMatrix<STATE> * TPZStructMatrixTBB::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *stiff = Create();
    
    int64_t cols = MAX(1, rhs.Cols());
    rhs.Redim(fEquationFilter.NEqExpand(),cols);
    
    Assemble(*stiff,rhs,guiInterface);
    
#ifdef LOG4CXX2
    if(loggerel->isDebugEnabled())
    {
        std::stringstream sout;
        stiff->Print("Stiffness matrix",sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
#endif
    return stiff;
    
}

void TPZStructMatrixTBB::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    this->fFlowGraph->ExecuteGraph(&rhs, &mat);
#else
    std::cout << "To use the tbb flow graph assemble please compile the NeoPZ with USING_TBB." << std::endl;
    DebugStop();
#endif
}


void TPZStructMatrixTBB::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    this->fFlowGraph->ExecuteGraph(&rhs);
#else
    std::cout << "To use the tbb flow graph assemble please compile the NeoPZ with USING_TBB." << std::endl;
    DebugStop();
#endif
}



static bool CanAssemble(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        if (elContribute[connectlist[i]] >= 0){
            return false;
        }
    }
    return true;
}

static void AssembleColor(int el,TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        elContribute[connectlist[i]] = el;
    }
}

static int WhoBlockedMe(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &elSeqinv)
{
    int el = -1;
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int elBlocked = elContribute[connectlist[i]];
        if (elBlocked == -1) continue;
        int elBlockedIndex = elSeqinv[elBlocked];
        if (el == -1) el = elBlockedIndex;
        if (elBlockedIndex < el) el = elBlockedIndex;
    }
    return el;
}

static void RemoveEl(int el,TPZCompMesh *cmesh,TPZVec<int64_t> &elContribute,int64_t elSequence)
{
    TPZCompEl *cel = cmesh->ElementVec()[el];
    if(!cel) DebugStop();
    TPZStack<int64_t> connectlist;
    cel->BuildConnectList(connectlist);
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int conindex = connectlist[i];
        if (elContribute[conindex] != elSequence){
            DebugStop();
        }
        elContribute[conindex] = -1;
    }
}

static int MinPassIndex(TPZStack<int64_t> &connectlist,TPZVec<int64_t> &elContribute, TPZVec<int64_t> &passIndex)
{
    int minPassIndex = -1;
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int elcont = elContribute[connectlist[i]];
        int passindex = -1;
        if (elcont != -1){
            passindex = passIndex[elcont];
            if (minPassIndex == -1) minPassIndex = passindex;
        }
        if (minPassIndex < passindex) minPassIndex = passindex;
    }
    return minPassIndex;
}

#ifdef USING_TBB
void TPZStructMatrixTBB::TPZFlowGraph::ElementColoring()
{
    
    const int64_t nnodes = fCMesh->NConnects();
    const int64_t nel = fCMesh->ElementVec().NElements();
    
    TPZManVector<int64_t> elContribute(nnodes,-1), passIndex(nel,-1);
    
    fFirstElColor.Push(0);
    felSequenceColor.Resize(nel);
    felSequenceColor.Fill(-1);
    felSequenceColorInv.Resize(nel, -1);
    felSequenceColorInv.Fill(-1);
    fnextBlocked.Resize(nel);
    fnextBlocked.Fill(-1);
    int nelProcessed = 0;
    int currentEl = 0;
    int currentPassIndex = 0;
    while (nelProcessed < fElementOrder.NElements()){
        
        int64_t elindex = fElementOrder[currentEl];
        
        if(felSequenceColorInv[elindex] == -1)
        {
            TPZCompEl *cel = fCMesh->ElementVec()[elindex];
            
            
            if(!cel) continue;
            TPZStack<int64_t> connectlist;
            cel->BuildConnectList(connectlist);
//                 std::cout << "elcontribute " << elContribute << std::endl;
//                 std::cout << "connectlist " << connectlist << std::endl;
            int minPass = MinPassIndex(connectlist,elContribute,passIndex);
            if (minPass == -1){
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                felSequenceColor[nelProcessed] = elindex;
                felSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else if (minPass == currentPassIndex){
            }
            else if (minPass < currentPassIndex){
                while (!CanAssemble(connectlist,elContribute)){
                    const int el = WhoBlockedMe(connectlist,elContribute, felSequenceColorInv);
                    if (fnextBlocked[el] == -1) fnextBlocked[el] = nelProcessed;
                    int locindex = felSequenceColor[el];
                    RemoveEl(locindex,fCMesh,elContribute,locindex);
//                             std::cout << "elcontribute " << elContribute << std::endl;
                }
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                felSequenceColor[nelProcessed] = elindex;
                felSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else{
                DebugStop();
            }
        }
        currentEl++;
        if (currentEl == fElementOrder.NElements()){
            fFirstElColor.Push(nelProcessed);
            currentEl = 0;
            currentPassIndex++;
        }
    }
    
#ifdef PZDEBUG
    std::ofstream toto("../ColorMeshDebug.txt");
    toto << "elSequence\n" << fElementOrder << std::endl;
    toto << "elSequenceColor\n" << felSequenceColor << std::endl;
    toto << "elSequenceColorInv\n" << felSequenceColorInv << std::endl;
    toto << "elBlocked\n" <<  fnextBlocked << std::endl;
    toto << "elContribute\n" << elContribute << std::endl;
    toto << "fFirstElColor " << fFirstElColor << std::endl;
    toto << "passIndex\n" << passIndex << std::endl;
    toto.close();
#endif
}

TPZStructMatrixTBB::TPZFlowGraph::TPZFlowGraph(TPZStructMatrixTBB *strmat, bool onlyrhs)
: fCMesh(strmat->Mesh()), fStruct(strmat), fGlobMatrix(0), fGlobRhs(0), fOnlyRhs(onlyrhs)
{
    this->OrderElements();
    this->ElementColoring();
    if (fOnlyRhs) {
        this->CreateGraphRhs();
    }
    else
    {
        this->CreateGraph();
    }
}

TPZStructMatrixTBB::TPZFlowGraph::~TPZFlowGraph()
{
}

void TPZStructMatrixTBB::TPZFlowGraph::OrderElements()
{
    int numelconnected = 0;
    int nconnect = fCMesh->ConnectVec().NElements();
    int ic;
//    firstelconnect contains the first element index in the elconnect vector
    TPZVec<int> firstelconnect(nconnect+1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        numelconnected += fCMesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic+1] = firstelconnect[ic]+fCMesh->ConnectVec()[ic].NElConnected();
    }
//    cout << "numelconnected " << numelconnected << endl;
//    cout << "firstelconnect ";
//     for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
    TPZVec<int> elconnect(numelconnected,-1);
    int el;
    TPZCompEl *cel;
    for(el=0; el<fCMesh->ElementVec().NElements(); el++) {
        cel = fCMesh->ElementVec()[el];
        if(!cel) continue;
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int nc = connectlist.NElements();
        int ic;
        for(ic=0; ic<nc; ic++) {
            int cindex = connectlist[ic];
            elconnect[firstelconnect[cindex]] = el;
            firstelconnect[cindex]++;
        }
    }
//     for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
    firstelconnect[0] = 0;
    
    for(ic=0; ic<nconnect; ic++) {
        firstelconnect[ic+1] = firstelconnect[ic]+fCMesh->ConnectVec()[ic].NElConnected();
    }
    
    fElementOrder.Resize(fCMesh->ElementVec().NElements(),-1);
    fElementOrder.Fill(-1);
    TPZVec<int> nodeorder(fCMesh->ConnectVec().NElements(),-1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        int seqnum = fCMesh->ConnectVec()[ic].SequenceNumber();
        if(seqnum >= 0) nodeorder[seqnum] = ic;
    }
//     cout << "nodeorder ";
    /*  for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
     cout << endl;
     cout.flush();*/
    
    int seq;
    int elsequence = 0;
    TPZVec<int> elorderinv(fCMesh->ElementVec().NElements(),-1);
    for(seq=0; seq<nconnect; seq++) {
        ic = nodeorder[seq];
        if(ic == -1) continue;
        int firstind = firstelconnect[ic];
        int lastind = firstelconnect[ic+1];
        int ind;
        for(ind=firstind; ind<lastind; ind++) {
            el = elconnect[ind];
            if(el == -1) {
                continue;
            }
            if(elorderinv[el]==-1) elorderinv[el] = elsequence++;
        }
    }
    elsequence = 0;
    for(seq=0;seq<fCMesh->ElementVec().NElements();seq++) {
        if(elorderinv[seq] == -1) continue;
        fElementOrder[elorderinv[seq]] = seq;
    }
    
    for(seq=0;seq<fCMesh->ElementVec().NElements();seq++) {
        if(fElementOrder[seq]==-1) break;
    }
    
    fElementOrder.Resize(seq);
}

TPZStructMatrixTBB::TPZFlowGraph::TPZFlowGraph(TPZFlowGraph const &copy)
: fCMesh(copy.fStruct->Mesh()), fStruct(copy.fStruct), fGlobMatrix(0), fGlobRhs(0), fOnlyRhs(copy.fOnlyRhs), fFirstElColor(copy.fFirstElColor)
{
    this->fnextBlocked = copy.fnextBlocked;
    this->felSequenceColor = copy.felSequenceColor;
    this->felSequenceColorInv = copy.felSequenceColorInv;
    this->fElementOrder = copy.fElementOrder;
    if (fOnlyRhs) {
        this->CreateGraphRhs();
    }
    else
    {
        this->CreateGraph();
    }
}

void TPZStructMatrixTBB::TPZFlowGraph::CreateGraphRhs()
{
    // create nodes for successive sum of rhs
    int64_t numcolors = fFirstElColor.size()-1;
    fNodeDest.resize(numcolors);
    fNodeDest.Fill(-1);
    int twoexp = 0;
    TPZStack<std::pair<int64_t,int64_t> > sumcolors;
    TPZStack<int64_t> numreceive;
    std::set<int64_t> received;
    int index = 0;
    while (numcolors>1) {
        int64_t num = 0;
        twoexp++;
        int64_t halfnumcolors = (numcolors / 2) + numcolors%2;
        for (int64_t i=halfnumcolors; i< numcolors; i++) {
            int64_t nr=0;
            if (received.find(i-halfnumcolors) == received.end()) {
                nr++;
                received.insert(i-halfnumcolors);
                fNodeDest[i-halfnumcolors] = sumcolors.size();
            }
            if (received.find(i) == received.end()) {
                nr++;
                received.insert(i);
                fNodeDest[i] = sumcolors.size();
            }
            std::pair<int64_t, int64_t> temp(i-halfnumcolors,i);
            sumcolors.Push(temp);
            numreceive.Push(nr);
        }
        numcolors = halfnumcolors;
    }
    numcolors = fFirstElColor.size()-1;
    fNodes.resize(sumcolors.size());
    // create the nodes and the links between them
    // the first nodes that mentions an index receives the rhs
    for (int64_t i=0; i<sumcolors.size(); i++) {
        TSumTwoColors block(sumcolors[i].first,sumcolors[i].second,&fRhsFat);
        fNodes[i] = new tbb::flow::continue_node<tbb::flow::continue_msg>(fGraph,numreceive[i],block);
    }
    for (int64_t i = sumcolors.size()-1; i>0; i--) {
        int rhsindex = sumcolors[i].first;
        for (int64_t j=i-1; j>=0; j++) {
            if (sumcolors[j].first == rhsindex || sumcolors[j].second == rhsindex) {
                tbb::flow::make_edge(*fNodes[i], *fNodes[j]);
                break;
            }
        }
    }
    int64_t neq = fCMesh->NEquations();
    fRhsFat.Redim(neq, numcolors);
    
}



void TPZStructMatrixTBB::TPZFlowGraph::CreateGraph()
{
    int64_t nelem = fCMesh->NElements();
    int64_t nconnects = fCMesh->NConnects();
    int64_t numberOfElements=felSequenceColor.NElements();
    this->felSequenceColor=felSequenceColor;
    
    TPZVec<int64_t> elementloaded(nconnects,-1);
    
    fNodes.resize(numberOfElements);
    fElMatPointers.Resize(numberOfElements);
    for (int64_t iel=0; iel<numberOfElements; iel++) {
        TPZAssembleTask body(iel,this);
        fNodes[iel] = new tbb::flow::continue_node<tbb::flow::continue_msg>(fGraph,1,body);
    }
    
    
    for (int64_t graphindex = 0; graphindex<numberOfElements; graphindex++) {
        int64_t el = felSequenceColor[graphindex];
        TPZCompEl *cel = fCMesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZStack<int64_t> connects;
        cel->BuildConnectList(connects);
        std::set<int64_t> fromwhere;
        for (int ic=0; ic<connects.size(); ic++) {
            int64_t c = connects[ic];
            if (elementloaded[c] != -1) {
                int64_t elorig = elementloaded[c];
//                in order to compute only once
                if (fromwhere.find(elorig) == fromwhere.end()) {
#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "Adding edge from " << elorig << " to " << graphindex;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    
                    make_edge(*fNodes[elorig], *fNodes[graphindex]);
                }
                fromwhere.insert(elorig);
            }
        }
        
        for (int ic=0; ic<connects.size(); ic++) {
            int64_t c = connects[ic];
            elementloaded[c] = graphindex;
        }
    }
    
}


void TPZStructMatrixTBB::TPZFlowGraph::ExecuteGraph(TPZFMatrix<STATE> *rhs, TPZMatrix<STATE> *matrix)
{
    if (fOnlyRhs && matrix != 0) {
        DebugStop();
    }
    this->fGlobMatrix = matrix;
    this->fGlobRhs = rhs;
    
    
    TPZCalcTask calcTasks(this);
    parallel_for(tbb::blocked_range<int64_t>(0, felSequenceColor.size()), calcTasks );
    
    fGraph.wait_for_all();
}

void TPZStructMatrixTBB::TPZFlowGraph::ExecuteGraph(TPZFMatrix<STATE> *rhs)
{
    if (!fOnlyRhs) {
        DebugStop();
    }
    this->fGlobRhs = rhs;
    this->fGlobMatrix = 0;
    int64_t numcolors = fFirstElColor.size()-1;
    fGlobRhs->Redim(this->fGlobRhs->Rows(), numcolors);
    
    TAssembleOneColor onecolor(this);
    parallel_for(tbb::blocked_range<int64_t>(0, fFirstElColor.size()-1), onecolor );
    
    fGraph.wait_for_all();
    int64_t nr = fRhsFat.Rows();
    for (int64_t r=0; r<nr; r++) {
        (*fGlobRhs)(r,0) = fRhsFat(r,0);
    }
}

void TPZStructMatrixTBB::TPZFlowGraph::TAssembleOneColor::operator()(const tbb::blocked_range<int64_t> &range) const
{
    for(int64_t color = range.begin(); color != range.end(); ++color)
    {
        TComputeElementRange elrange(fFlowGraph,color);
        int64_t firstel = fFlowGraph->fFirstElColor[color];
        int64_t lastel = fFlowGraph->fFirstElColor[color+1];
        parallel_for(tbb::blocked_range<int64_t>(firstel,lastel),elrange);
        // trigger the sum node
        int64_t node = fFlowGraph->fNodeDest[color];
        (fFlowGraph->fNodes)[node]->try_put(tbb::flow::continue_msg());
        
    }
}

void TPZStructMatrixTBB::TPZFlowGraph::TComputeElementRange::operator()(const tbb::blocked_range<int64_t> &range) const
{
    TPZCompMesh *cmesh = fFlowGraph->fCMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = fFlowGraph->fGuiInterface;
    TPZElementMatrix ef(cmesh, TPZElementMatrix::EF);
    
    for (int64_t iel=range.begin(); iel != range.end(); iel++)
    {
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "Computing element " << iel;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef LOG4CXX
        std::stringstream sout;
        sout << "Element " << iel << " elapsed time ";
        TPZTimer timeforel(sout.str());
        timeforel.start();
#endif
        
        int element = fFlowGraph->felSequenceColor[iel];
        
        if (element >= 0){
            
            TPZCompEl *el = cmesh->ElementVec()[element];
            
            if(!el) return;
            
            el->CalcResidual(ef);
            
            if(!el->HasDependency()) {
                
                ef.ComputeDestinationIndices();
                fFlowGraph->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            } else
            {
                ef.ApplyConstraints();
                ef.ComputeDestinationIndices();
                fFlowGraph->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            }
            
        }
        
        int64_t nrows = fFlowGraph->fGlobRhs->Rows();
        TPZFMatrix<STATE> locrhs(nrows,1,&fFlowGraph->fRhsFat(0, fColor),nrows);
        
        if(!ef.HasDependency()) {
            
            locrhs.AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
        } else
        {
            locrhs.AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
        }
        
    } // outsided if
    
    
}


tbb::flow::continue_msg TPZStructMatrixTBB::TPZFlowGraph::TPZAssembleTask::operator()(const tbb::flow::continue_msg &msg)
{
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " Element " << this->fIel;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZElementMatrix *Ek = fElMat->fEk;
    TPZElementMatrix *Ef = fElMat->fEf;
#ifdef PZDEBUG
    if (!Ef) {
        DebugStop();
    }
#endif
    
    if(fOrigin->fGlobMatrix) {
//        assemble the matrix
        if(!Ek->HasDependency()) {
            fOrigin->fGlobMatrix->AddKel(Ek->fMat,Ek->fSourceIndex,Ek->fDestinationIndex);
            fOrigin->fGlobRhs->AddFel(Ef->fMat,Ek->fSourceIndex,Ek->fDestinationIndex);
        } else {
            fOrigin->fGlobMatrix->AddKel(Ek->fConstrMat,Ek->fSourceIndex,Ek->fDestinationIndex);
            fOrigin->fGlobRhs->AddFel(Ef->fConstrMat,Ek->fSourceIndex,Ek->fDestinationIndex);
        }
    } else {
        if(!Ef->HasDependency()) {
            fOrigin->fGlobRhs->AddFel(Ef->fMat,Ef->fSourceIndex,Ef->fDestinationIndex);
        } else {
            fOrigin->fGlobRhs->AddFel(Ef->fConstrMat,Ef->fSourceIndex,Ef->fDestinationIndex);
        }
    }
    if (Ek) {
        delete Ek;
        fElMat->fEk = 0;
    }
    delete Ef;
    fElMat->fEf = 0;
    
    return tbb::flow::continue_msg();
}

void TPZStructMatrixTBB::TPZFlowGraph::TPZCalcTask::operator()(const tbb::blocked_range<int64_t>& range) const
{
    TPZCompMesh *cMesh = fFlowGraph->fCMesh;
    
    TPZVec<int64_t> &elSequenceColor = fFlowGraph->felSequenceColor;
    
    for(int iel = range.begin(); iel != range.end(); ++iel) {
        
        int element = elSequenceColor[iel];
        
        if (element < 0)
        {
            DebugStop();
        }
        
        TPZAssembleTask Assemble = tbb::flow::copy_body<TPZAssembleTask,tbb::flow::continue_node<tbb::flow::continue_msg> >(*fFlowGraph->fNodes[element]);
        if (Assemble.fIel != element) {
            DebugStop();
        }
        if (fFlowGraph->fGlobMatrix)
        {
            (Assemble.fElMat->fEk) = new TPZElementMatrix(cMesh,TPZElementMatrix::EK);
        }
        Assemble.fElMat->fEf = new TPZElementMatrix(cMesh,TPZElementMatrix::EF);
        
        TPZElementMatrix *ek = (Assemble.fElMat->fEk);
        TPZElementMatrix *ef = (Assemble.fElMat->fEf);
        
        TPZCompEl *el = cMesh->ElementVec()[element];
        
        if(!el) continue;
        
        if (fFlowGraph->fGlobMatrix)
            el->CalcStiff(*ek,*ef);
        else
            el->CalcResidual(*ef);
        
        if(!el->HasDependency()) {
            
            if (fFlowGraph->fGlobMatrix) {
                ek->ComputeDestinationIndices();
                fFlowGraph->fStruct->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
            } else {
                ef->ComputeDestinationIndices();
                fFlowGraph->fStruct->FilterEquations(ef->fSourceIndex,ef->fDestinationIndex);
            }
            
        } else {
//            the element has dependent nodes
            if (fFlowGraph->fGlobMatrix) {
                ek->ApplyConstraints();
                ef->ApplyConstraints();
                ek->ComputeDestinationIndices();
                fFlowGraph->fStruct->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
            } else {
                ef->ApplyConstraints();
                ef->ComputeDestinationIndices();
                fFlowGraph->fStruct->FilterEquations(ef->fSourceIndex,ef->fDestinationIndex);
            }
        }
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " Element " << element;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        (fFlowGraph->fNodes)[element]->try_put(tbb::flow::continue_msg());
        
    }
} // operator()

#endif

int TPZStructMatrixTBB::ClassId() const{
    return Hash("TPZStructMatrixTBB") ^ TPZStructMatrixBase::ClassId() << 1;
}

void TPZStructMatrixTBB::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf, context);
#ifdef USING_TBB
    fFlowGraph = dynamic_cast<TPZFlowGraph *>(TPZPersistenceManager::GetInstance(&buf));
#endif
}

void TPZStructMatrixTBB::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);
    DebugStop();
#ifdef USING_TBB
//    TPZPersistenceManager::WritePointer(fFlowGraph, &buf);
#endif
}

//void TPZStructMatrixTBB::TPZFlowNode::operator()(tbb::flow::continue_msg) const
//{
//   TPZCompMesh *cmesh = myGraph->fStruct->Mesh();
//   TPZAutoPointer<TPZGuiInterface> guiInterface = myGraph->fGuiInterface;
//   TPZElementMatrix ek(cmesh, TPZElementMatrix::EK);
//   TPZElementMatrix ef(cmesh, TPZElementMatrix::EF);
//#ifdef LOG4CXX
//   if (logger->isDebugEnabled()) {
//       std::stringstream sout;
//       sout << "Computing element " << iel;
//       LOGPZ_DEBUG(logger, sout.str())
//   }
//#endif
//#ifdef LOG4CXX
//   std::stringstream sout;
//   sout << "Element " << iel << " elapsed time ";
//   TPZTimer timeforel(sout.str());
//   timeforel.start();
//#endif
//
//   int element = myGraph->felSequenceColor[iel];
//
//   if (element >= 0){
//
//       TPZCompEl *el = cmesh->ElementVec()[element];
//
//       if(!el) return;
//
//       if (myGraph->fGlobMatrix)
//           el->CalcStiff(ek,ef);
//       else
//           el->CalcResidual(ef);
//
//       if(!el->HasDependency()) {
//
//           if (myGraph->fGlobMatrix) {
//               ek.ComputeDestinationIndices();
//               myGraph->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
//           } else {
//               ef.ComputeDestinationIndices();
//               myGraph->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
//           }
//
//       } else {
//           // the element has dependent nodes
//           if (myGraph->fGlobMatrix) {
//               ek.ApplyConstraints();
//               ef.ApplyConstraints();
//               ek.ComputeDestinationIndices();
//               myGraph->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
//           } else {
//               ef.ApplyConstraints();
//               ef.ComputeDestinationIndices();
//               myGraph->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
//           }
//
//       }
//
//
//       if(myGraph->fGlobMatrix) {
//           // assemble the matrix
//           if(!ek.HasDependency()) {
//               myGraph->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
//               myGraph->fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
//           } else {
//               myGraph->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
//               myGraph->fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
//           }
//       } else {
//           if(!ef.HasDependency()) {
//               myGraph->fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
//           } else {
//               myGraph->fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
//           }
//       }
//
//   } // outsided if
//
//#ifdef LOG4CXX
//   timeforel.stop();
//   if (logger->isDebugEnabled())
//   {
//       std::stringstream sout;
//       sout << timeforel.processName() <<  timeforel;
//       LOGPZ_DEBUG(logger, sout.str())
//   }
//#endif
//
//}

template class TPZRestoreClass<TPZStructMatrixTBB>;
#endif
