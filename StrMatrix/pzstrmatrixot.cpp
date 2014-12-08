/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixOT methods.
 */

#include "pzstrmatrixot.h"

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
#include "run_stats_table.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixOT"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph_ot("-ass_graph_ot", "Run statistics table for the graph creation and coloring TPZStructMatrixOT.");


TPZStructMatrixOT::TPZStructMatrixOT(TPZCompMesh *mesh) : fMesh(mesh), fEquationFilter(mesh->NEquations()) {
    fMesh = mesh;
    this->SetNumThreads(0);
    stat_ass_graph_ot.start();
    TPZManVector<int> ElementOrder;
    TPZStructMatrixOT::OrderElement(this->Mesh(), ElementOrder);
    TPZStructMatrixOT::ElementColoring(this->Mesh(), ElementOrder, felSequenceColor, fnextBlocked, ColorsToProcess);
    stat_ass_graph_ot.stop();
}

TPZStructMatrixOT::TPZStructMatrixOT(TPZAutoPointer<TPZCompMesh> cmesh) : fCompMesh(cmesh), fEquationFilter(cmesh->NEquations()) {
    fMesh = cmesh.operator->();
    this->SetNumThreads(0);
    
    stat_ass_graph_ot.start();
    TPZManVector<int> ElementOrder;
    TPZStructMatrixOT::OrderElement(this->Mesh(), ElementOrder);
    TPZStructMatrixOT::ElementColoring(this->Mesh(), ElementOrder, felSequenceColor, fnextBlocked, ColorsToProcess);
    stat_ass_graph_ot.stop();
}

TPZStructMatrixOT::TPZStructMatrixOT(const TPZStructMatrixOT &copy) : fMesh(copy.fMesh), fEquationFilter(copy.fEquationFilter)
{
    if (copy.fCompMesh) {
        fCompMesh = copy.fCompMesh;
    }
    fMaterialIds = copy.fMaterialIds;
    fNumThreads = copy.fNumThreads;
    felSequenceColor = copy.felSequenceColor;
    fnextBlocked = copy.fnextBlocked;
    ColorsToProcess = copy.ColorsToProcess;
}

TPZMatrix<STATE> *TPZStructMatrixOT::Create() {
    cout << "TPZStructMatrixOT::Create should never be called\n";
    return 0;
}

TPZStructMatrixOT *TPZStructMatrixOT::Clone() {
    cout << "TPZStructMatrixOT::Clone should never be called\n";
    return 0;
}

RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixOT::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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

void TPZStructMatrixOT::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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



void TPZStructMatrixOT::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ){
    
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

void TPZStructMatrixOT::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
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
void TPZStructMatrixOT::FilterEquations(TPZVec<long> &origindex, TPZVec<long> &destindex) const
{
    //destindex = origindex;
    fEquationFilter.Filter(origindex, destindex);
    
}

TPZMatrix<STATE> * TPZStructMatrixOT::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *stiff = Create();
    
    //long neq = stiff->Rows();
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
void TPZStructMatrixOT::SetMaterialIds(const std::set<int> &materialids)
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


void TPZStructMatrixOT::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    ThreadData threaddata(this,mat,rhs,fMaterialIds,guiInterface);
    
    threaddata.fnextBlocked=&fnextBlocked;
    threaddata.felSequenceColor=&felSequenceColor;
    
    const int numthreads = this->fNumThreads;
    TPZVec<pthread_t> allthreads(numthreads);
    int itr;
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return;
        }
    }
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWork,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
    
#ifdef LOG4CXX
    if(loggerCheck->isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        mat.Print("Matriz de Rigidez: ",sout,EMathematicaInput);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
}


void TPZStructMatrixOT::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    ThreadData threaddata(this,rhs,fMaterialIds,guiInterface);
    
    threaddata.fnextBlocked=&fnextBlocked;
    threaddata.felSequenceColor=&felSequenceColor;
    
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return;
        }
    }
    
#ifdef USING_TBB
    int begin=0;
    
    for(int color=0; color<ColorsToProcess.NElements(); color++) {
        tbb::task_group tg;
        for (int t=begin; t<ColorsToProcess[color]; t++) {
            tg.run(WorkResidualTBB(t, &threaddata));
        }
        tg.wait();
        begin=ColorsToProcess[color];
    }
#endif
    
}



TPZStructMatrixOT::ThreadData::ThreadData(TPZStructMatrixOT *strmat, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZStructMatrixOT::ThreadData::ThreadData()");
    pthread_cond_init(&fCondition, NULL);
    /*	sem_t *sem_open( ... );
     int sem_close(sem_t *sem);
     int sem_unlink(const char *name);
     */
    /*
     #ifdef MACOSX
     std::stringstream sout;
     static int counter = 0;
     sout << "AssemblySem" << counter++;
     fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
     if(fAssembly == SEM_FAILED)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     DebugStop();
     }
     #else
     int sem_result = sem_init(&fAssembly,0,0);
     if(sem_result != 0)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     }
     #endif
     */
}

TPZStructMatrixOT::ThreadData::ThreadData(TPZStructMatrixOT *strmat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat),fGuiInterface(guiInterface), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZStructMatrixOT::ThreadData::ThreadData()");
    pthread_cond_init(&fCondition, NULL);
    /*	sem_t *sem_open( ... );
     int sem_close(sem_t *sem);
     int sem_unlink(const char *name);
     */
    /*
     #ifdef MACOSX
     std::stringstream sout;
     static int counter = 0;
     sout << "AssemblySem" << counter++;
     fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
     if(fAssembly == SEM_FAILED)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     DebugStop();
     }
     #else
     int sem_result = sem_init(&fAssembly,0,0);
     if(sem_result != 0)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     }
     #endif
     */
}

TPZStructMatrixOT::ThreadData::~ThreadData()
{
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZStructMatrixOT::ThreadData::~ThreadData()");
    pthread_cond_destroy(&fCondition);
    /*
     #ifdef MACOSX
     sem_close(fAssembly);
     #else
     sem_destroy(&fAssembly);
     #endif
     */
}

void *TPZStructMatrixOT::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    int nel = cmesh->NElements();
    bool hasWork = true;
    int iel = 0;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    while(hasWork)
    {
        tht::EnterCriticalSection(data->fAccessElement);
        
        int localiel = data->fNextElement;
        if(data->fNextElement < nel){
            if (!data->felBlocked.size() || data->felBlocked.begin()->first > localiel){
                iel = (*data->felSequenceColor)[localiel];
                
                if(localiel==-1) DebugStop();
                
                int elBl = (*data->fnextBlocked)[localiel];
                if (elBl >= 0){
                    data->felBlocked[elBl]++;
                    hasWork = true;
                }
                
                data->fNextElement++;
            }
            else if (data->felBlocked.size() || data->felBlocked.begin()->first <= localiel){
                data->fSleeping=true;
                while(data->fSleeping){
                    pthread_cond_wait(&data->fCondition, &data->fAccessElement);
                }
                iel = -1;
                hasWork = true;
            }
            else{
                DebugStop();
            }
        }
        else{
            hasWork = false;
            iel = -1;
        }
        
        tht::LeaveCriticalSection(data->fAccessElement);
        
        
        if (iel >= 0){
            TPZCompEl *el = cmesh->ElementVec()[iel];
            el->CalcStiff(ek,ef);
            if(guiInterface) if(guiInterface->AmIKilled()){
                break;
            }
            if(!el->HasDependency()) {
                
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef LOG4CXX
                if(loggerel->isDebugEnabled())
                {
                    std::stringstream sout;
                    ek.fMat.Print("Element stiffness matrix",sout);
                    ef.fMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            } else {
                // the element has dependent nodes
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef LOG4CXX
                if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
                {
                    std::stringstream sout;
                    el->Reference()->Print(sout);
                    el->Print(sout);
                    ek.Print(sout);
                    ef.Print(sout);
                    LOGPZ_DEBUG(loggerel2,sout.str())
                }
#endif
            }
            
            
            if(data->fGlobMatrix){
                // Assemble the matrix
                if(!ek.HasDependency())
                {
                    data->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                    data->fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                else
                {
                    data->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                    data->fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                
            }
            
            tht::EnterCriticalSection( data->fAccessElement );
            
            int elBl = (*data->fnextBlocked)[localiel];
            if (elBl >= 0 && data->felBlocked.find(elBl) != data->felBlocked.end()){
                data->felBlocked[elBl]--;
                if (data->felBlocked[elBl] == 0){
                    data->felBlocked.erase(elBl);
                    if(data->fSleeping) {
                        data->fSleeping=false;
                    }
                    pthread_cond_broadcast(&data->fCondition);
                }
            }
            else if(elBl >= 0){
                DebugStop();
            }
            
            tht::LeaveCriticalSection( data->fAccessElement );
        }
    }
    return 0;
}

void *TPZStructMatrixOT::ThreadData::ThreadWorkResidual(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    int nel = cmesh->NElements();
    bool hasWork = true;
    int iel = 0;
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    while(hasWork)
    {
        tht::EnterCriticalSection(data->fAccessElement);
        int localiel = data->fNextElement;
        if(data->fNextElement < nel){
            if (!data->felBlocked.size() || data->felBlocked.begin()->first > localiel){
                iel = (*data->felSequenceColor)[localiel];
                
                if(localiel==-1) DebugStop();
                
                int elBl = (*data->fnextBlocked)[localiel];
                if (elBl >= 0){
                    data->felBlocked[elBl]++;
                    hasWork = true;
                }
                
                data->fNextElement++;
            }
            else if (data->felBlocked.size() || data->felBlocked.begin()->first <= localiel){
                data->fSleeping=true;
                while(data->fSleeping){
                    pthread_cond_wait(&data->fCondition, &data->fAccessElement);
                }
                iel = -1;
                hasWork = true;
            }
            else{
                DebugStop();
            }
        }
        else{
            hasWork = false;
            iel = -1;
        }
        
        tht::LeaveCriticalSection(data->fAccessElement);
        if (iel >= 0){
            TPZCompEl *el = cmesh->ElementVec()[iel];
            
            el->CalcResidual(ef);
            
            if(guiInterface) if(guiInterface->AmIKilled()){
                break;
            }
            if(!el->HasDependency()) {
                ef.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                
#ifdef LOG4CXX
                if(loggerel->isDebugEnabled())
                {
                    std::stringstream sout;
                    ef.fMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            } else {
                // the element has dependent nodes
                ef.ApplyConstraints();
                ef.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
#ifdef LOG4CXX
                if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
                {
                    std::stringstream sout;
                    el->Reference()->Print(sout);
                    el->Print(sout);
                    LOGPZ_DEBUG(loggerel2,sout.str())
                }
#endif
            }
            
            if(data->fGlobRhs){
                // Assemble the matrix
                if(!ef.HasDependency())
                {
                    data->fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
                }
                else
                {
                    data->fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
                }
            }
            tht::EnterCriticalSection( data->fAccessElement );
            int elBl = (*data->fnextBlocked)[localiel];
            if (elBl >= 0 && data->felBlocked.find(elBl) != data->felBlocked.end()){
                data->felBlocked[elBl]--;
                if (data->felBlocked[elBl] == 0){
                    data->felBlocked.erase(elBl);
                    if(data->fSleeping) {
                        data->fSleeping=false;
                    }
                    pthread_cond_broadcast(&data->fCondition);
                }
            }
            else if(elBl >= 0){
                DebugStop();
            }
            tht::LeaveCriticalSection(data->fAccessElement);
        }
    }
    
    return 0;
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

void TPZStructMatrixOT::ElementColoring(TPZCompMesh *cmesh, TPZVec<int> &elSequence, TPZVec<int> &elSequenceColor,
                                        TPZVec<int> &elBlocked, TPZVec<int> &elColors)
{
    
    const int nnodes = cmesh->NConnects();
    const int nel = cmesh->ElementVec().NElements();
    
    if (!nel) return;
    
    elColors.Resize(nel, -1);
    
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
            elColors[currentPassIndex]=nelProcessed;
            currentPassIndex++;
        }
    }
    
    elColors[currentPassIndex]=elColors[currentPassIndex-1]+1;
    
    elColors.Resize(currentPassIndex+1);
    
    
    
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

void TPZStructMatrixOT::OrderElement(TPZCompMesh *cmesh, TPZVec<int> &ElementOrder)
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

TPZStructMatrixOT::WorkResidualTBB::WorkResidualTBB(int elem, ThreadData *data)
:fElem(elem), data(data)
{
    
    
    
}
void TPZStructMatrixOT::WorkResidualTBB::operator()()
{
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    
    int element = (*data->felSequenceColor)[fElem];
    
    if (element >= 0) {
        
        TPZCompEl *el = cmesh->ElementVec()[element];
        
        if (data->fGlobMatrix) {
            el->CalcStiff(ek,ef);
        } else {
            el->CalcResidual(ef);
        }
        
        if(!el->HasDependency()) {
            ef.ComputeDestinationIndices();
            data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
            
            if (data->fGlobMatrix) {
                ek.ApplyConstraints();
                ek.ComputeDestinationIndices();
            }
        }
        
        if(!ef.HasDependency()) {
            data->fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
            if (data->fGlobMatrix) {
                data->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            }
        }
        else {
            data->fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
            if (data->fGlobMatrix) {
                data->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            }
        }
    }
}

#endif