#include "TPZStructMatrixTBBFlowUtils.h"
#include "pzstrmatrixflowtbb.h"
#include "TPZStructMatrix.h"
#include "pzcmesh.h"
#include "TPZTimer.h"
#include "TPZElementMatrixT.h"
#include "pzlog.h"

#ifdef USING_TBB

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixTBBFlow");
#endif

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

template<class TVar>
void TPZFlowGraph<TVar>::ElementColoring()
{
    
    const int nnodes = fMesh->NConnects();
    const int nel = fMesh->ElementVec().NElements();
    
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
            TPZCompEl *cel = fMesh->ElementVec()[elindex];
            
            
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
                    RemoveEl(locindex,fMesh,elContribute,locindex);
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

template<class TVar>
TPZFlowGraph<TVar>::TPZFlowGraph(TPZStructMatrixTBBFlow<TVar> *strmat, bool computeRhs)
    : fStartNode(fGraph), fStruct(strmat),
      fGlobMatrix(0), fGlobRhs(0), fComputeRhs(computeRhs)
{
//    if(
//        !dynamic_cast<TPZMatrix<TVar>*>(&stiff_base)||
//        !dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)
//       ){
//        PZError<<__PRETTY_FUNCTION__;
//        PZError<<" Incompatible type. Aborting...\n";
//        DebugStop();
//    }
    auto *myself = dynamic_cast<TPZStructMatrix*>(strmat);
    fMesh = myself->Mesh();

    this->OrderElements();
    this->ElementColoring();
    this->CreateGraph();
}

template<class TVar>
TPZFlowGraph<TVar>::~TPZFlowGraph()
{
    for (int k = 0; k < fNodes.size(); ++k) {
        delete fNodes[k];
    }
    
}

template<class TVar>
void TPZFlowGraph<TVar>::ExecuteGraph(TPZBaseMatrix *rhs, TPZBaseMatrix *matrix)
{
    
    this->fGlobMatrix = dynamic_cast<TPZMatrix<TVar>*>(matrix);
    this->fGlobRhs = dynamic_cast<TPZFMatrix<TVar>*>(rhs);
    if(!fGlobMatrix || !fGlobRhs)
    {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" Incompatible type. Aborting...\n";
        DebugStop();
    }
    this->fStartNode.try_put(tbb::flow::continue_msg());
    this->fGraph.wait_for_all();
    
}

template<class TVar>
TPZFlowGraph<TVar>::TPZFlowGraph(const TPZFlowGraph<TVar> &copy)
: fStartNode(fGraph), fStruct(copy.fStruct), fGlobMatrix(0), fGlobRhs(0)
{
    auto *myself = dynamic_cast<TPZStructMatrix*>(fStruct);
    fMesh = myself->Mesh();

    this->fnextBlocked = copy.fnextBlocked;
    this->felSequenceColor = copy.felSequenceColor;
    this->felSequenceColorInv = copy.felSequenceColorInv;
    this->fElementOrder = copy.fElementOrder;
    this->CreateGraph();
}

template<class TVar>
void TPZFlowNode<TVar>::operator()(tbb::flow::continue_msg) const
{
    auto *mystruct = dynamic_cast<TPZStructMatrix*>(myGraph->fStruct);
    TPZCompMesh *cmesh = mystruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = myGraph->fGuiInterface;
    TPZElementMatrixT<TVar> ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrixT<TVar> ef(cmesh,TPZElementMatrix::EF);
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
                mystruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
            } else {
                ef.ComputeDestinationIndices();
                mystruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            }
            
        } else {
            // the element has dependent nodes
            if (myGraph->fGlobMatrix) {
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                mystruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
            } else {
                ef.ApplyConstraints();
                ef.ComputeDestinationIndices();
                mystruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            }
            
        }
        
        
        if(myGraph->fGlobMatrix) {
            // assemble the matrix
            if(!ek.HasDependency()) {
                myGraph->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                if(fComputeRhs)
                    myGraph->fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            } else {
                myGraph->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                if(fComputeRhs)
                    myGraph->fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            }
        } else {
            if(fComputeRhs && !ef.HasDependency()) {
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

template<class TVar>
void TPZFlowGraph<TVar>::CreateGraph()
{
    int64_t nelem = fMesh->NElements();
    int64_t nconnects = fMesh->NConnects();
    int64_t numberOfElements=felSequenceColor.NElements();
    
    // each graphnode represents an element that can be computed and assembled
    fNodes.resize(felSequenceColor.NElements());
    for (int64_t i=0; i<felSequenceColor.NElements(); i++) {
        fNodes[i]= new tbb::flow::continue_node<tbb::flow::continue_msg>(fGraph, TPZFlowNode(this, i, fComputeRhs));
    }
    TPZVec<int64_t> elementloaded(nconnects,-1);
    
    for (int64_t graphindex = 0; graphindex<numberOfElements; graphindex++) {
        int64_t el = felSequenceColor[graphindex];
        TPZCompEl *cel = fMesh->Element(el);
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

template<class TVar>
void TPZFlowGraph<TVar>::OrderElements()
{
    int numelconnected = 0;
    int nconnect = fMesh->ConnectVec().NElements();
    int ic;
    //firstelconnect contains the first element index in the elconnect vector
    TPZVec<int> firstelconnect(nconnect+1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        numelconnected += fMesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic+1] = firstelconnect[ic]+fMesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "numelconnected " << numelconnected << endl;
    //cout << "firstelconnect ";
    //  for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
    TPZVec<int> elconnect(numelconnected,-1);
    int el;
    TPZCompEl *cel;
    for(el=0; el<fMesh->ElementVec().NElements(); el++) {
        cel = fMesh->ElementVec()[el];
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
        firstelconnect[ic+1] = firstelconnect[ic]+fMesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "elconnect\n";
    //  int no;
    //  for(no=0; no< fMesh->ConnectVec().NElements(); no++) {
    //cout << "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
    //       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
    //cout << endl;
    //  }
    
    fElementOrder.Resize(fMesh->ElementVec().NElements(),-1);
    fElementOrder.Fill(-1);
    TPZVec<int> nodeorder(fMesh->ConnectVec().NElements(),-1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        int seqnum = fMesh->ConnectVec()[ic].SequenceNumber();
        if(seqnum >= 0) nodeorder[seqnum] = ic;
    }
    //  cout << "nodeorder ";
    /*  for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
     cout << endl;
     cout.flush();*/
    int seq;
    int elsequence = 0;
    TPZVec<int> elorderinv(fMesh->ElementVec().NElements(),-1);
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
    for(seq=0;seq<fMesh->ElementVec().NElements();seq++) {
        if(elorderinv[seq] == -1) continue;
        fElementOrder[elorderinv[seq]] = seq;
    }
    
    for(seq=0;seq<fMesh->ElementVec().NElements();seq++) {
        if(fElementOrder[seq]==-1) break;
    }
    
    fElementOrder.Resize(seq);
}
        
template class TPZFlowGraph<STATE>;
template class TPZFlowNode<STATE>;
template class TPZFlowGraph<CSTATE>;
template class TPZFlowNode<CSTATE>;
#endif
