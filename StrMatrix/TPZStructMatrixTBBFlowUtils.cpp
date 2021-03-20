#include "TPZStructMatrixTBBFlowUtils.h"
#include "pzstrmatrixflowtbb.h"
#include "pzcmesh.h"
#include "TPZTimer.h"
#include "pzlog.h"

#ifdef USING_TBB

static TPZLogger logger("pz.strmatrix.TPZStructMatrixTBBFlow");

static int MinPassIndex(TPZStack<int64_t> &connectlist,TPZVec<int> &elContribute, TPZVec<int> &passIndex)
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

static void AssembleColor(int el,TPZStack<int64_t> &connectlist, TPZVec<int> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        elContribute[connectlist[i]] = el;
    }
}

static bool CanAssemble(TPZStack<int64_t> &connectlist, TPZVec<int> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        if (elContribute[connectlist[i]] >= 0){
            return false;
        }
    }
    return true;
}

static int WhoBlockedMe(TPZStack<int64_t> &connectlist, TPZVec<int> &elContribute, TPZVec<int> &elSeqinv)
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

static void RemoveEl(int el,TPZCompMesh *cmesh,TPZVec<int> &elContribute,int elSequence)
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

void TPZFlowGraph::ElementColoring()
{
    
    const int nnodes = cmesh->NConnects();
    const int nel = cmesh->ElementVec().NElements();
    
    TPZManVector<int> elContribute(nnodes,-1), passIndex(nel,-1);
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
        
        int elindex = fElementOrder[currentEl];
        
        if(felSequenceColorInv[elindex] == -1)
        {
            TPZCompEl *cel = cmesh->ElementVec()[elindex];
            
            
            if(!cel) continue;
            TPZStack<int64_t> connectlist;
            cel->BuildConnectList(connectlist);
            //      std::cout << "elcontribute " << elContribute << std::endl;
            //      std::cout << "connectlist " << connectlist << std::endl;
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
                    RemoveEl(locindex,cmesh,elContribute,locindex);
                    //          std::cout << "elcontribute " << elContribute << std::endl;
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
            currentEl = 0;
            currentPassIndex++;
        }
    }
    
    //std::cout << "sequence: " << elSequence << std::endl;
    //std::cout << "color: " << elSequenceColorInv << std::endl;
    
    
    //    exit(101);
#ifdef PZDEBUG
    std::ofstream toto("../ColorMeshDebug.txt");
    toto << "elSequence\n" << fElementOrder << std::endl;
    toto << "elSequenceColor\n" << felSequenceColor << std::endl;
    toto << "elSequenceColorInv\n" << felSequenceColorInv << std::endl;
    toto << "elBlocked\n" <<  fnextBlocked << std::endl;
    toto << "elContribute\n" << elContribute << std::endl;
    toto << "passIndex\n" << passIndex << std::endl;
    toto.close();
#endif
}

TPZFlowGraph::TPZFlowGraph(TPZStructMatrixTBBFlow *strmat)
: cmesh(strmat->Mesh()), fStartNode(fGraph), fStruct(strmat), fGlobMatrix(0), fGlobRhs(0)
{
    this->OrderElements();
    this->ElementColoring();
    this->CreateGraph();
}

TPZFlowGraph::~TPZFlowGraph()
{
    for (int k = 0; k < fNodes.size(); ++k) {
        delete fNodes[k];
    }
    
}

void TPZFlowGraph::ExecuteGraph(TPZFMatrix<STATE> *rhs, TPZMatrix<STATE> *matrix)
{
    
    this->fGlobMatrix = matrix;
    this->fGlobRhs = rhs;
    this->fStartNode.try_put(tbb::flow::continue_msg());
    this->fGraph.wait_for_all();
    
}

TPZFlowGraph::TPZFlowGraph(TPZFlowGraph const &copy)
: cmesh(copy.fStruct->Mesh()), fStartNode(fGraph), fStruct(copy.fStruct), fGlobMatrix(0), fGlobRhs(0)
{
    this->fnextBlocked = copy.fnextBlocked;
    this->felSequenceColor = copy.felSequenceColor;
    this->felSequenceColorInv = copy.felSequenceColorInv;
    this->fElementOrder = copy.fElementOrder;
    this->CreateGraph();
}

void TPZFlowNode::operator()(tbb::flow::continue_msg) const
{
    TPZCompMesh *cmesh = myGraph->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = myGraph->fGuiInterface;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Computing element " << iel;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef PZ_LOG
    std::stringstream sout;
    sout << "Element " << iel << " elapsed time ";
    TPZTimer timeforel(sout.str());
    timeforel.start();
#endif
    
    int element = myGraph->felSequenceColor[iel];
    
    if (element >= 0){
        
        TPZCompEl *el = cmesh->ElementVec()[element];
        
        if(!el) return;
        
        if (myGraph->fGlobMatrix)
            el->CalcStiff(ek,ef);
        else
            el->CalcResidual(ef);
        
        if(!el->HasDependency()) {
            
            if (myGraph->fGlobMatrix) {
                ek.ComputeDestinationIndices();
                myGraph->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
            } else {
                ef.ComputeDestinationIndices();
                myGraph->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            }
            
        } else {
            // the element has dependent nodes
            if (myGraph->fGlobMatrix) {
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                myGraph->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
            } else {
                ef.ApplyConstraints();
                ef.ComputeDestinationIndices();
                myGraph->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            }
            
        }
        
        
        if(myGraph->fGlobMatrix) {
            // assemble the matrix
            if(!ek.HasDependency()) {
                myGraph->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                myGraph->fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            } else {
                myGraph->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                myGraph->fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            }
        } else {
            if(!ef.HasDependency()) {
                myGraph->fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
            } else {
                myGraph->fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
            }
        }
        
    } // outsided if
    
#ifdef PZ_LOG
    timeforel.stop();
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << timeforel.processName() <<  timeforel;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
}

void TPZFlowGraph::CreateGraph()
{
    int64_t nelem = cmesh->NElements();
    int64_t nconnects = cmesh->NConnects();
    int64_t numberOfElements=felSequenceColor.NElements();
    
    // each graphnode represents an element that can be computed and assembled
    fNodes.resize(felSequenceColor.NElements());
    for (int64_t i=0; i<felSequenceColor.NElements(); i++) {
        fNodes[i]= new tbb::flow::continue_node<tbb::flow::continue_msg>(fGraph, TPZFlowNode(this, i));
    }
    TPZVec<int64_t> elementloaded(nconnects,-1);
    
    for (int64_t graphindex = 0; graphindex<numberOfElements; graphindex++) {
        int64_t el = felSequenceColor[graphindex];
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZStack<int64_t> connects;
        cel->BuildConnectList(connects);
        int ngraphs = 0;
        std::set<int64_t> fromwhere;
        for (int ic=0; ic<connects.size(); ic++) {
            int64_t c = connects[ic];
            if (elementloaded[c] != -1) {
                int64_t elorig = elementloaded[c];
                // in order to compute only once
                if (fromwhere.find(elorig) == fromwhere.end()) {
#ifdef PZ_LOG
                    if (logger.isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "Adding edge from " << elorig << " to " << graphindex;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    make_edge(*fNodes[elorig], *fNodes[graphindex]);
                }
                fromwhere.insert(elorig);
                ngraphs++;
            }
        }
        if (ngraphs == 0) {
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Setting start element " << graphindex;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
            make_edge(fStartNode, *fNodes[graphindex]);
        }
        for (int ic=0; ic<connects.size(); ic++) {
            int64_t c = connects[ic];
            elementloaded[c] = graphindex;
        }
    }
    
}

void TPZFlowGraph::OrderElements()
{
    int numelconnected = 0;
    int nconnect = cmesh->ConnectVec().NElements();
    int ic;
    //firstelconnect contains the first element index in the elconnect vector
    TPZVec<int> firstelconnect(nconnect+1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        numelconnected += cmesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic+1] = firstelconnect[ic]+cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "numelconnected " << numelconnected << endl;
    //cout << "firstelconnect ";
    //  for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
    TPZVec<int> elconnect(numelconnected,-1);
    int el;
    TPZCompEl *cel;
    for(el=0; el<cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
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
    //  for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        firstelconnect[ic+1] = firstelconnect[ic]+cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "elconnect\n";
    //  int no;
    //  for(no=0; no< fMesh->ConnectVec().NElements(); no++) {
    //cout << "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
    //       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
    //cout << endl;
    //  }
    
    fElementOrder.Resize(cmesh->ElementVec().NElements(),-1);
    fElementOrder.Fill(-1);
    TPZVec<int> nodeorder(cmesh->ConnectVec().NElements(),-1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        int seqnum = cmesh->ConnectVec()[ic].SequenceNumber();
        if(seqnum >= 0) nodeorder[seqnum] = ic;
    }
    //  cout << "nodeorder ";
    /*  for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
     cout << endl;
     cout.flush();*/
    int seq;
    int elsequence = 0;
    TPZVec<int> elorderinv(cmesh->ElementVec().NElements(),-1);
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
    //  cout << "elorderinv ";
    //  for(seq=0;seq<fMesh->ElementVec().NElements();seq++) cout << elorderinv[seq] << ' ';
    //  cout << endl;
    elsequence = 0;
    for(seq=0;seq<cmesh->ElementVec().NElements();seq++) {
        if(elorderinv[seq] == -1) continue;
        fElementOrder[elorderinv[seq]] = seq;
    }
    
    for(seq=0;seq<cmesh->ElementVec().NElements();seq++) {
        if(fElementOrder[seq]==-1) break;
    }
    
    fElementOrder.Resize(seq);
}
#endif