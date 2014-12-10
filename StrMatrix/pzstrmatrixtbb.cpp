/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixTBB methods.
 */

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
#include "pzmaterial.h"

using namespace std;

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

#ifdef USING_TBB
using namespace tbb::flow;
#endif

#include "run_stats_table.h"

RunStatsTable stat_ass_graph_tbb("-ass_graph_tbb", "Run statistics table for the graph creation, coloring and tbb::flow::graph TPZStructMatrixTBB.");


TPZStructMatrixTBB::TPZStructMatrixTBB(TPZCompMesh *mesh) : fMesh(mesh), fEquationFilter(mesh->NEquations()) {
    fMesh = mesh;
    this->SetNumThreads(0);
#ifdef USING_TBB
    stat_ass_graph_tbb.start();
    TPZManVector<int> ElementOrder;
    TPZStructMatrixTBB::OrderElement(this->Mesh(), ElementOrder);
    TPZStructMatrixTBB::ElementColoring(this->Mesh(), ElementOrder, felSequenceColor, fnextBlocked);
    stat_ass_graph_tbb.stop();
#endif
}

TPZStructMatrixTBB::TPZStructMatrixTBB(TPZAutoPointer<TPZCompMesh> cmesh) : fCompMesh(cmesh), fEquationFilter(cmesh->NEquations()) {
    fMesh = cmesh.operator->();
    this->SetNumThreads(0);
#ifdef USING_TBB
    stat_ass_graph_tbb.start();
    TPZManVector<int> ElementOrder;
    TPZStructMatrixTBB::OrderElement(this->Mesh(), ElementOrder);
    TPZStructMatrixTBB::ElementColoring(this->Mesh(), ElementOrder, felSequenceColor, fnextBlocked);
    stat_ass_graph_tbb.stop();
#endif
}

TPZStructMatrixTBB::TPZStructMatrixTBB(const TPZStructMatrixTBB &copy) : fMesh(copy.fMesh), fEquationFilter(copy.fEquationFilter)
{
    if (copy.fCompMesh) {
        fCompMesh = copy.fCompMesh;
    }
    fMaterialIds = copy.fMaterialIds;
    fNumThreads = copy.fNumThreads;
    felSequenceColor = copy.felSequenceColor;
    fnextBlocked = copy.fnextBlocked;

}

TPZStructMatrixTBB::~TPZStructMatrixTBB()
{
#ifdef USING_TBB
    if (fAssembleThreadGraph)
        delete fAssembleThreadGraph;
#endif
}
TPZMatrix<STATE> *TPZStructMatrixTBB::Create() {
    cout << "TPZStructMatrixTBB::Create should never be called\n";
    return 0;
}

TPZStructMatrixTBB *TPZStructMatrixTBB::Clone() {
    cout << "TPZStructMatrixTBB::Clone should never be called\n";
    return 0;
}

RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixTBB::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        long neqcondense = fEquationFilter.NActiveEquations();
#ifdef DEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense,rhs.Cols(),0.);
        if(this->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhsloc,guiInterface);
        }
        else{
            this->Serial_Assemble(stiffness,rhsloc,guiInterface);
        }
        
        fEquationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        if(this->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhs,guiInterface);
        }
        else{
            this->Serial_Assemble(stiffness,rhs,guiInterface);
        }
    }
    ass_stiff.stop();
}

void TPZStructMatrixTBB::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_rhs.start();
    if(fEquationFilter.IsActive())
    {
        long neqcondense = fEquationFilter.NActiveEquations();
        long neqexpand = fEquationFilter.NEqExpand();
        if(rhs.Rows() != neqexpand || Norm(rhs) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense,1,0.);
        if(this->fNumThreads)
        {
            this->MultiThread_Assemble(rhsloc,guiInterface);
        }
        else
        {
            this->Serial_Assemble(rhsloc,guiInterface);
        }
        fEquationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        if(this->fNumThreads){
            this->MultiThread_Assemble(rhs,guiInterface);
        }
        else{
            this->Serial_Assemble(rhs,guiInterface);
        }
    }
    ass_rhs.stop();
}



void TPZStructMatrixTBB::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ){
    
    if(!fMesh){
        LOGPZ_ERROR(logger,"Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef LOG4CXX
    if(dynamic_cast<TPZSubCompMesh * >(fMesh))
    {
        std::stringstream sout;
        sout << "AllEig = {};";
        LOGPZ_DEBUG(loggerelmat,sout.str())
        
    }
#endif
#ifdef DEBUG
    if (rhs.Rows() != fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif
    
    long iel;
    long nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK),ef(fMesh, TPZElementMatrix::EF);
#ifdef LOG4CXX
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
    long count = 0;
    for(iel=0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        int matidsize = fMaterialIds.size();
        if(matidsize){
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat)
            {
                if (!submesh) {
                    continue;
                }
                else if(submesh->NeedsComputing(fMaterialIds) == false) continue;
            }
            else
            {
                int matid = mat->Id();
                if (this->ShouldCompute(matid) == false) continue;
            }
        }
        
        count++;
        if(!(count%1000))
        {
            std::cout << '*';
            std::cout.flush();
        }
        if(!(count%20000))
        {
            std::cout << "\n";
        }
        calcstiff.start();
        
        el->CalcStiff(ek,ef);
        
        if(guiInterface) if(guiInterface->AmIKilled()){
            return;
        }
        
#ifdef LOG4CXX
        if(dynamic_cast<TPZSubCompMesh * >(fMesh))
        {
            std::stringstream objname;
            objname << "Element" << iel;
            std::string name = objname.str();
            objname << " = ";
            std::stringstream sout;
            ek.fMat.Print(objname.str().c_str(),sout,EMathematicaInput);
            sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";
            
            LOGPZ_DEBUG(loggerelmat,sout.str())
            /*		  if(iel == 133)
             {
             std::stringstream sout2;
             el->Reference()->Print(sout2);
             el->Print(sout2);
             LOGPZ_DEBUG(logger,sout2.str())
             }
             */
        }
        
#endif
        
#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if(!result)
        {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck,sout.str())
        }
        
#endif
        
        calcstiff.stop();
        assemble.start();
        
        if(!el->HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //			TPZSFMatrix<STATE> test(stiffness);
            //			TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //			stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //			rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //			test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			test -= stiffness;
            //			test.Print("matriz de rigidez diference",std::cout);
            //			test2.Print("matriz de rigidez interface",std::cout);
            
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled())
            {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if(gel)
                {
                    TPZManVector<REAL> center(gel->Dimension()),xcenter(3,0.);
                    gel->CenterPoint(gel->NSides()-1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                }
                else {
                    sout << "Stiffness for computational element without associated geometric element\n";
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled() && ! dynamic_cast<TPZSubCompMesh *>(fMesh))
            {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                TPZManVector<REAL> center(gel->Dimension()),xcenter(3,0.);
                gel->CenterPoint(gel->NSides()-1, center);
                gel->X(center, xcenter);
                sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        }
        
        assemble.stop();
    }//fim for iel
    if(count > 20) std::cout << std::endl;
    
#ifdef LOG4CXX
    if(loggerCheck->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        stiffness.Print("Matriz de Rigidez: ",sout,EMathematicaInput);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
    
#endif
    
}

void TPZStructMatrixTBB::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
    long iel;
    long nelem = fMesh->NElements();
    
    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");
    
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
    for(iel=0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        
        TPZMaterial * mat = el->Material();
        if (!mat) continue;
        int matid = mat->Id();
        if (this->ShouldCompute(matid) == false) continue;
        
        TPZElementMatrix ef(fMesh, TPZElementMatrix::EF);
        
        calcresidual.start();
        
        el->CalcResidual(ef);
        
        calcresidual.stop();
        
        assemble.start();
        
        if(!el->HasDependency()) {
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
        }
        
        assemble.stop();
        
    }//fim for iel
#ifdef LOG4CXX
    {
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << calcresidual.processName() << " " << calcresidual << std::endl;
            sout << assemble.processName() << " " << assemble;
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
    }
#endif
    //std::cout << std::endl;
}

/// filter out the equations which are out of the range
void TPZStructMatrixTBB::FilterEquations(TPZVec<long> &origindex, TPZVec<long> &destindex) const
{
    //destindex = origindex;
    fEquationFilter.Filter(origindex, destindex);
    
}

TPZMatrix<STATE> * TPZStructMatrixTBB::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *stiff = Create();

    long cols = MAX(1, rhs.Cols());
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

/// Set the set of material ids which will be considered when assembling the system
void TPZStructMatrixTBB::SetMaterialIds(const std::set<int> &materialids)
{
    fMaterialIds = materialids;
#ifdef LOG4CXX
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
    long iel;
    TPZAdmChunkVector<TPZCompEl*> &elvec = fMesh->ElementVec();
    long nel = elvec.NElements();
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
        TPZAutoPointer<TPZStructMatrix> str = anal->StructMatrix();
        if(!str)
        {
            LOGPZ_WARN(logger,"SetMaterialIds called for substructure without structural matrix")
            continue;
        }
        str->SetMaterialIds(materialids);
    }
}


void TPZStructMatrixTBB::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    
    TPZGraphThreadData AssembleThreadGraph(fnextBlocked, felSequenceColor);
    
    AssembleThreadGraph.fStruct=this;
    AssembleThreadGraph.fGlobMatrix=&mat;
    AssembleThreadGraph.fGlobRhs=&rhs;
    AssembleThreadGraph.fStart.try_put(continue_msg());
    AssembleThreadGraph.fAssembleGraph.wait_for_all();
    
#endif
}


void TPZStructMatrixTBB::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    TPZGraphThreadData AssembleThreadGraph(fnextBlocked, felSequenceColor);
    
    AssembleThreadGraph.fStruct=this;
    AssembleThreadGraph.fGlobMatrix=0;
    AssembleThreadGraph.fGlobRhs=&rhs;
    AssembleThreadGraph.fStart.try_put(continue_msg());
    AssembleThreadGraph.fAssembleGraph.wait_for_all();
    
#endif
}



static bool CanAssemble(TPZStack<long> &connectlist, TPZVec<int> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        if (elContribute[connectlist[i]] >= 0){
            return false;
        }
    }
    return true;
}

static void AssembleColor(int el,TPZStack<long> &connectlist, TPZVec<int> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        elContribute[connectlist[i]] = el;
    }
}

static int WhoBlockedMe(TPZStack<long> &connectlist, TPZVec<int> &elContribute, TPZVec<int> &elSeqinv)
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
    TPZStack<long> connectlist;
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

static int MinPassIndex(TPZStack<long> &connectlist,TPZVec<int> &elContribute, TPZVec<int> &passIndex)
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

void TPZStructMatrixTBB::ElementColoring(TPZCompMesh *cmesh, TPZVec<int> &elSequence, TPZVec<int> &elSequenceColor,
                                         TPZVec<int> &elBlocked)
{
    
    const int nnodes = cmesh->NConnects();
    const int nel = cmesh->ElementVec().NElements();
    
    TPZManVector<int> elContribute(nnodes,-1), passIndex(nel,-1), elSequenceColorInv(nel,-1);
    elSequenceColor.Resize(nel);
    elSequenceColor.Fill(-1);
    elBlocked.Resize(nel);
    elBlocked.Fill(-1);
    int nelProcessed = 0;
    int currentEl = 0;
    int currentPassIndex = 0;
    while (nelProcessed < elSequence.NElements()){
        
        int elindex = elSequence[currentEl];
        
        if(elSequenceColorInv[elindex] == -1)
        {
            TPZCompEl *cel = cmesh->ElementVec()[elindex];
            
            
            if(!cel) continue;
            TPZStack<long> connectlist;
            cel->BuildConnectList(connectlist);
            //      std::cout << "elcontribute " << elContribute << std::endl;
            //      std::cout << "connectlist " << connectlist << std::endl;
            int minPass = MinPassIndex(connectlist,elContribute,passIndex);
            if (minPass == -1){
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                elSequenceColor[nelProcessed] = elindex;
                elSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else if (minPass == currentPassIndex){
            }
            else if (minPass < currentPassIndex){
                while (!CanAssemble(connectlist,elContribute)){
                    const int el = WhoBlockedMe(connectlist,elContribute, elSequenceColorInv);
                    if (elBlocked[el] == -1) elBlocked[el] = nelProcessed;
                    int locindex = elSequenceColor[el];
                    RemoveEl(locindex,cmesh,elContribute,locindex);
                    //          std::cout << "elcontribute " << elContribute << std::endl;
                }
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                elSequenceColor[nelProcessed] = elindex;
                elSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else{
                DebugStop();
            }
        }
        currentEl++;
        if (currentEl == elSequence.NElements()){
            currentEl = 0;
            currentPassIndex++;
        }
    }
    
    //std::cout << "sequence: " << elSequence << std::endl;
    //std::cout << "color: " << elSequenceColorInv << std::endl;
    
    
    //    exit(101);
    /*
     std::ofstream toto("c:\\Temp\\output\\ColorMeshDebug.txt");
     toto << "elSequence\n" << elSequence << std::endl;
     toto << "elSequenceColor\n" << elSequenceColor << std::endl;
     toto << "elSequenceColorInv\n" << elSequenceColorInv << std::endl;
     toto << "elBlocked\n" << elBlocked << std::endl;
     toto << "elContribute\n" << elContribute << std::endl;
     toto << "passIndex\n" << passIndex << std::endl;
     toto.close();
     */
}

void TPZStructMatrixTBB::OrderElement(TPZCompMesh *cmesh, TPZVec<int> &ElementOrder)
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
        TPZStack<long> connectlist;
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
    
    ElementOrder.Resize(cmesh->ElementVec().NElements(),-1);
    ElementOrder.Fill(-1);
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
        ElementOrder[elorderinv[seq]] = seq;
    }
    
    for(seq=0;seq<cmesh->ElementVec().NElements();seq++) {
        if(ElementOrder[seq]==-1) break;
    }
    
    ElementOrder.Resize(seq);
}

#ifdef USING_TBB

void TPZStructMatrixTBB::TPZGraphThreadNode::operator()(tbb::flow::continue_msg) const{
    
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    
    int element = data->felSequenceColor[iel];
    
            if (element >= 0){
            
            TPZCompEl *el = cmesh->ElementVec()[element];
                
            if(!el) return;
            
            if (data->fGlobMatrix)
                el->CalcStiff(ek,ef);
            else
                el->CalcResidual(ef);
            
            if(!el->HasDependency()) {
                
                if (data->fGlobMatrix) {
                    ek.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                } else {
                    ef.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                }
                
            } else {
                // the element has dependent nodes
                if (data->fGlobMatrix) {
                    ek.ApplyConstraints();
                    ef.ApplyConstraints();
                    ek.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                } else {
                    ef.ApplyConstraints();
                    ef.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                }
                
            }
            
            
            if(data->fGlobMatrix) {
                // assemble the matrix
                if(!ek.HasDependency()) {
                    data->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                    data->fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                } else {
                    data->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                    data->fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
            } else {
                if(!ef.HasDependency()) {
                    data->fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
                } else {
                    data->fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
                }
            }
                
        } // outsided if
        
}

// destructor
TPZStructMatrixTBB::TPZGraphThreadData::~TPZGraphThreadData() {
    for (int i=0; i<fGraphNodes.size(); i++) {
        delete fGraphNodes[i];
    }
}

TPZStructMatrixTBB::TPZGraphThreadData::TPZGraphThreadData(TPZVec<int> fnextBlocked, TPZVec<int> felSequenceColor)
: fStart(fAssembleGraph)
{
    this->felSequenceColor=felSequenceColor;
    fGraphNodes.resize(felSequenceColor.NElements());
    for (int i=0; i<felSequenceColor.NElements(); i++) {
        fGraphNodes[i]= new continue_node<continue_msg>(fAssembleGraph, TPZGraphThreadNode(this, i));
    }
    // initial edges to independent nodes
    for (int i=0; i<felSequenceColor.NElements(); i++) {
        make_edge(fStart, *fGraphNodes[i]);
    }
    // edges from sucessors to nextBlockedElements
    for (int i=0; i<felSequenceColor.NElements(); i++) {
        int next=fnextBlocked[i];
        if (next>0) make_edge(*fGraphNodes[i], *fGraphNodes[next]);
    }
 
}

#endif