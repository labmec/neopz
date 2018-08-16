/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixGC methods.
 */

#include "pzstrmatrixgc.h"

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

using namespace std;

#include "pzlog.h"

#include "pz_pthread.h"
#include "run_stats_table.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixGC"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph("-ass_graph", "Run statistics table for the graph creation and coloring TPZStructMatrixGC.");


TPZStructMatrixGC::TPZStructMatrixGC(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {
    stat_ass_graph.start();
    TPZManVector<int64_t> ElementOrder;
    TPZStructMatrixGC::OrderElement(this->Mesh(), ElementOrder);
    TPZStructMatrixGC::ElementColoring(this->Mesh(), ElementOrder, felSequenceColor, fnextBlocked);
    stat_ass_graph.stop();
}

TPZStructMatrixGC::TPZStructMatrixGC(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
    stat_ass_graph.start();
    TPZManVector<int64_t> ElementOrder;
    TPZStructMatrixGC::OrderElement(this->Mesh(), ElementOrder);
    TPZStructMatrixGC::ElementColoring(this->Mesh(), ElementOrder, felSequenceColor, fnextBlocked);
    stat_ass_graph.stop();
}

TPZStructMatrixGC::TPZStructMatrixGC(const TPZStructMatrixGC &copy) : TPZStructMatrixBase(copy){
    felSequenceColor = copy.felSequenceColor;
    fnextBlocked = copy.fnextBlocked;
}

TPZMatrix<STATE> *TPZStructMatrixGC::Create() {
    cout << "TPZStructMatrixGC::Create should never be called\n";
    return 0;
}

TPZStructMatrixGC *TPZStructMatrixGC::Clone() {
    cout << "TPZStructMatrixGC::Clone should never be called\n";
    return 0;
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");


void TPZStructMatrixGC::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
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

void TPZStructMatrixGC::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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



void TPZStructMatrixGC::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ){
    
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
#ifdef PZDEBUG
    if (rhs.Rows() != fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif
    
    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK),ef(fMesh, TPZElementMatrix::EF);
#ifdef LOG4CXX
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
    int64_t count = 0;
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

void TPZStructMatrixGC::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
    int64_t iel;
    int64_t nelem = fMesh->NElements();
    
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

TPZMatrix<STATE> * TPZStructMatrixGC::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *stiff = Create();
    
    //int64_t neq = stiff->Rows();
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

void TPZStructMatrixGC::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
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


void TPZStructMatrixGC::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    ThreadData threaddata(this,rhs,fMaterialIds,guiInterface);
    
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
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWorkResidual,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
}



TPZStructMatrixGC::ThreadData::ThreadData(TPZStructMatrixGC *strmat, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZStructMatrixGC::ThreadData::ThreadData()");
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

TPZStructMatrixGC::ThreadData::ThreadData(TPZStructMatrixGC *strmat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat),fGuiInterface(guiInterface), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZStructMatrixGC::ThreadData::ThreadData()");
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

TPZStructMatrixGC::ThreadData::~ThreadData()
{
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZStructMatrixGC::ThreadData::~ThreadData()");
    pthread_cond_destroy(&fCondition);
    /*
     #ifdef MACOSX
     sem_close(fAssembly);
     #else
     sem_destroy(&fAssembly);
     #endif
     */
}

void *TPZStructMatrixGC::ThreadData::ThreadWork(void *datavoid)
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
        // nextelement is a protected value
        int64_t localiel = data->fNextElement;
        int64_t blockedel = -1;
        if(data->felBlocked.size()) blockedel = data->felBlocked.begin()->first;
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Checking out " << localiel << " next blocked element " << blockedel;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        if(localiel < nel)
        {
            // The elblocked data structure indicates the elements blocked by the elements being processed (why is this needed?)
            if (!data->felBlocked.size() || data->felBlocked.begin()->first > localiel)
            {
                
                if(localiel==-1) DebugStop();
                // this is the next element that will be computed
                iel = (*data->felSequenceColor)[localiel];
                
                // identify the element which will be blocked by iel
                int elBl = (*data->fnextBlocked)[localiel];
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element can be computed iel = " << localiel << " element blocked " << elBl;
                    LOGPZ_DEBUG(logger, sout.str())
                }
                
#endif
                // update the datastructure with the highest element that can be processed while iel hasn't been computed yet
                if (elBl >= 0){
                    data->felBlocked[elBl]++;
                    hasWork = true;
                }
                // the next thread will monitor a higher numbered element
                data->fNextElement++;
            }
            else if (data->felBlocked.size() || data->felBlocked.begin()->first <= localiel){
#ifdef LOG4CXX
                if (logger->isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "Going to sleep cannot do " << localiel << " because of " << data->felBlocked.begin()->first << " has to be computed first";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                // localiel cannot be processed yet
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
        
#ifdef LOG4CXX
        std::stringstream sout;
        sout << "Element " << localiel << " elapsed time ";
        TPZTimer timeforel(sout.str());
        timeforel.start();
#endif
        
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
            
#ifdef LOG4CXX
            timeforel.stop();
            if (logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << timeforel.processName() <<  timeforel;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
            tht::EnterCriticalSection( data->fAccessElement );
            
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Computed and Assembled " << localiel;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            // clean up the data structure
            int elBl = (*data->fnextBlocked)[localiel];
            if (elBl >= 0 && data->felBlocked.find(elBl) != data->felBlocked.end())
            {
                data->felBlocked[elBl]--;
                int dataelbl = data->felBlocked[elBl];
                if (dataelbl == 0){
#ifdef LOG4CXX
                    if(logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Element " << elBl << " released";
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    data->felBlocked.erase(elBl);
                    if(data->fSleeping) {
                        data->fSleeping=false;
                    
#ifdef LOG4CXX
                        if(logger->isDebugEnabled()){
                            LOGPZ_DEBUG(logger, "Waking up everybody")
                        }
#endif
                        // wake up everybody
                        pthread_cond_broadcast(&data->fCondition);
                    }
                }
                else if (dataelbl < 0)
                {
                    DebugStop();
                }
                else
                {
#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "Not freeing the blocked element " << elBl;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
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

void *TPZStructMatrixGC::ThreadData::ThreadWorkResidual(void *datavoid)
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

static bool CanAssemble(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute) {
    for (int64_t i = 0; i < connectlist.NElements(); i++) {
        if (elContribute[connectlist[i]] >= 0) {
            return false;
        }
    }
    return true;
}

static void AssembleColor(int64_t el, TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute) {
    for (int64_t i = 0; i < connectlist.NElements(); i++) {
        elContribute[connectlist[i]] = el;
    }
}

static int64_t WhoBlockedMe(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &elSeqinv) {
    int64_t el = -1;
    for (int64_t i = 0; i < connectlist.NElements(); i++) {
        int64_t elBlocked = elContribute[connectlist[i]];
        if (elBlocked == -1) continue;
        int64_t elBlockedIndex = elSeqinv[elBlocked];
        if (el == -1) el = elBlockedIndex;
        if (elBlockedIndex < el) el = elBlockedIndex;
    }
    return el;
}

static void RemoveEl(int64_t el, TPZCompMesh *cmesh, TPZVec<int64_t> &elContribute, int64_t elSequence) {
    TPZCompEl *cel = cmesh->ElementVec()[el];
    if (!cel) DebugStop();
    TPZStack<int64_t> connectlist;
    cel->BuildConnectList(connectlist);
    for (int64_t i = 0; i < connectlist.NElements(); i++) {
        int64_t conindex = connectlist[i];
        if (elContribute[conindex] != elSequence) {
            DebugStop();
        }
        elContribute[conindex] = -1;
    }
}

static int MinPassIndex(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &passIndex) {
    int minPassIndex = -1;
    for (int64_t i = 0; i < connectlist.NElements(); i++) {
        int64_t elcont = elContribute[connectlist[i]];
        int64_t passindex = -1;
        if (elcont != -1) {
            passindex = passIndex[elcont];
            if (minPassIndex == -1) minPassIndex = passindex;
        }
        if (minPassIndex < passindex) minPassIndex = passindex;
    }
    return minPassIndex;
}

void TPZStructMatrixGC::ElementColoring(TPZCompMesh *cmesh, TPZVec<int64_t> &elSequence, TPZVec<int64_t> &elSequenceColor,
        TPZVec<int64_t> &elBlocked) {

    const int64_t n_connects = cmesh->NConnects();
    const int64_t nel = cmesh->NElements();

    TPZManVector<int64_t> elContribute(n_connects, -1); // given a connect index, tells the index of the element which will assemble it
    TPZManVector<int64_t> elSequenceColorInv(nel, -1); // given an element index, tells its color
    TPZManVector<int64_t> passIndex(nel, -1); // the number of the pass in which the element was colored
    elSequenceColor.Resize(nel);
    elSequenceColor.Fill(-1);
    elBlocked.Resize(nel);
    elBlocked.Fill(-1);
    int64_t nelProcessed = 0; // Total number of colored elements so far
    int64_t currentPassIndex = 0; // Index of this pass (a certain number of passes through the list of elements is needed for full coloring)
    while (nelProcessed < elSequence.NElements()) {
        for (auto elindex : elSequence) {
            if (elSequenceColorInv[elindex] == -1) {
                TPZCompEl *cel = cmesh->Element(elindex);

                if (!cel) continue;
                TPZStack<int64_t> connectlist;
                cel->BuildConnectList(connectlist);
                //      std::cout << "elcontribute " << elContribute << std::endl;
                //      std::cout << "connectlist " << connectlist << std::endl;
                int minPass = MinPassIndex(connectlist, elContribute, passIndex);
                // None of the connects in this element has been associated with a colored element
                if (minPass == -1) {
                    passIndex[elindex] = currentPassIndex;
                    AssembleColor(elindex, connectlist, elContribute);
                    elSequenceColor[nelProcessed] = elindex;
                    elSequenceColorInv[elindex] = nelProcessed;
                    nelProcessed++;
                } else if (minPass == currentPassIndex) {
                } else if (minPass < currentPassIndex) {
                    while (!CanAssemble(connectlist, elContribute)) {
                        const int64_t el = WhoBlockedMe(connectlist, elContribute, elSequenceColorInv);
                        if (elBlocked[el] == -1) elBlocked[el] = nelProcessed;
                        int64_t locindex = elSequenceColor[el];
                        RemoveEl(locindex, cmesh, elContribute, locindex);
                        //          std::cout << "elcontribute " << elContribute << std::endl;
                    }
                    passIndex[elindex] = currentPassIndex;
                    AssembleColor(elindex, connectlist, elContribute);
                    elSequenceColor[nelProcessed] = elindex;
                    elSequenceColorInv[elindex] = nelProcessed;
                    nelProcessed++;
                } else {
                    DebugStop();
                }
            }
        }
        currentPassIndex++;
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

void TPZStructMatrixGC::OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder) {
    int64_t numelconnected = 0;
    int64_t nconnect = cmesh->NConnects();
    //firstelconnect contains the first element index in the elconnect vector
    TPZVec<int64_t> firstelconnect(nconnect + 1);
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        numelconnected += cmesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic + 1] = firstelconnect[ic] + cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "numelconnected " << numelconnected << endl;
    //cout << "firstelconnect ";
    //  for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
    TPZVec<int64_t> elconnect(numelconnected, -1);
    TPZCompEl *cel;
    for (int64_t el = 0; el < cmesh->NElements(); el++) {
        cel = cmesh->Element(el);
        if (!cel) continue;
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int64_t nc = connectlist.NElements();
        for (int64_t ic = 0; ic < nc; ic++) {
            int64_t cindex = connectlist[ic];
            elconnect[firstelconnect[cindex]] = el;
            firstelconnect[cindex]++;
        }
    }
    //  for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        firstelconnect[ic + 1] = firstelconnect[ic] + cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "elconnect\n";
    //  int no;
    //  for(no=0; no< fMesh->ConnectVec().NElements(); no++) {
    //cout << "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
    //       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
    //cout << endl;
    //  }

    ElementOrder.Resize(cmesh->ElementVec().NElements(), -1);
    ElementOrder.Fill(-1);
    TPZVec<int64_t> nodeorder(cmesh->ConnectVec().NElements(), -1);
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        int64_t seqnum = cmesh->ConnectVec()[ic].SequenceNumber();
        if (seqnum >= 0) nodeorder[seqnum] = ic;
    }
    //  cout << "nodeorder ";
    /*  for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
     cout << endl;
     cout.flush();*/
    int64_t elsequence = 0;
    TPZVec<int64_t> elorderinv(cmesh->ElementVec().NElements(), -1);
    for (int64_t seq = 0; seq < nconnect; seq++) {
        int64_t ic = nodeorder[seq];
        if (ic == -1) continue;
        int64_t firstind = firstelconnect[ic];
        int64_t lastind = firstelconnect[ic + 1];
        for (int64_t ind = firstind; ind < lastind; ind++) {
            int64_t el = elconnect[ind];
            if (el == -1) {
                continue;
            }
            if (elorderinv[el] == -1) elorderinv[el] = elsequence++;
        }
    }
    //  cout << "elorderinv ";
    //  for(seq=0;seq<fMesh->ElementVec().NElements();seq++) cout << elorderinv[seq] << ' ';
    //  cout << endl;
    elsequence = 0;
    for (int64_t seq = 0; seq < cmesh->ElementVec().NElements(); seq++) {
        if (elorderinv[seq] == -1) continue;
        ElementOrder[elorderinv[seq]] = seq;
    }
    
    int64_t seq;
    for (seq = 0; seq < cmesh->ElementVec().NElements(); seq++) {
        if (ElementOrder[seq] == -1) break;
    }
    ElementOrder.Resize(seq);
}

int TPZStructMatrixGC::ClassId() const{
    return Hash("TPZStructMatrixGC") ^ TPZStructMatrixBase::ClassId() << 1;
}

void TPZStructMatrixGC::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf, context);
    buf.Read(fnextBlocked);
    buf.Read(felSequenceColor);
}

void TPZStructMatrixGC::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);
    buf.Write(fnextBlocked);
    buf.Write(felSequenceColor);
}

template class TPZRestoreClass<TPZStructMatrixGC>;