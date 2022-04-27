/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixLCC methods.
 */

#include "pzstrmatrixLCC.h"

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
#include "TPZParallelUtils.h"
#include "pzgnode.h"
#include "TPZTimer.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#endif
#ifdef USING_OMP
#include "omp.h"
#endif
#include "pzcheckconsistency.h"
#include "TPZMaterial.h"

using namespace std;

#include "pzlog.h"

#include "run_stats_table.h"
#include <thread>

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixLCC");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerel2("pz.strmatrix.elementinterface");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph_LCC("-ass_graph_ot", "Run statistics table for the graph creation and coloring TPZStructMatrixLCC.");


TPZStructMatrixLCC::TPZStructMatrixLCC(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {

    //TPZStructMatrixLCC::OrderElements();
    
}

TPZStructMatrixLCC::TPZStructMatrixLCC(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
        
}

TPZStructMatrixLCC::TPZStructMatrixLCC(const TPZStructMatrixLCC &copy) : TPZStructMatrixBase(copy) {
    fElSequenceColor = copy.fElSequenceColor;
    fElBlocked = copy.fElBlocked;
    fEquationFilter = copy.fEquationFilter;
    fElementsComputed = copy.fElementsComputed;
    fElementCompleted = copy.fElementCompleted;
    fSomeoneIsSleeping = copy.fSomeoneIsSleeping;
    fShouldColor = copy.fShouldColor;
    fUsingTBB = copy.fUsingTBB;
}


TPZMatrix<STATE> *TPZStructMatrixLCC::Create() {
    cout << "TPZStructMatrixLCC::Create should never be called\n";
    return 0;
}

TPZStructMatrixLCC *TPZStructMatrixLCC::Clone() {
    cout << "TPZStructMatrixLCC::Clone should never be called\n";
    return 0;
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixLCC::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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

void TPZStructMatrixLCC::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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
#ifdef HUGEDEBUG
            TPZFMatrix<STATE> rhsserial(rhsloc);
            this->Serial_Assemble(rhsserial, guiInterface);
#endif
            this->MultiThread_Assemble(rhsloc,guiInterface);
#ifdef HUGEDEBUG
            rhsserial -= rhsloc;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
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
#ifdef HUGEDEBUG
            TPZFMatrix<STATE> rhsserial(rhs);
            this->Serial_Assemble(rhsserial, guiInterface);
#endif
            this->MultiThread_Assemble(rhs,guiInterface);
#ifdef HUGEDEBUG
            REAL normrhs = Norm(rhs);
            REAL normrhsserial = Norm(rhsserial);
            std::cout << "normrhs = " << normrhs << " normrhsserial " << normrhsserial << std::endl;
            rhsserial -= rhs;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
        }
        else{
            this->Serial_Assemble(rhs,guiInterface);
        }
    }
    ass_rhs.stop();
}



void TPZStructMatrixLCC::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ){
    int64_t nelem = fMesh->NElements();
                    
        for (int64_t iel = 0; iel < nelem; iel++)
        {
            TPZCompEl *el = fMesh->Element(iel);
            if (!el) continue;
            
            TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
            
            el->CalcStiff(ek, ef);
                        
            if(!el->HasDependency()) {
                ek.ComputeDestinationIndices();
                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
                //            TPZSFMatrix<STATE> test(stiffness);
                //            TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
                //            stiffness.Print("before assembly",std::cout,EMathematicaInput);
                stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                //            stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
                //            rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
                //            test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                //            test -= stiffness;
                //            test.Print("matriz de rigidez diference",std::cout);
                //            test2.Print("matriz de rigidez interface",std::cout);
                
    
            } else {
                // the element has dependent nodes
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
                stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
    
            }
        }
#ifdef PZDEBUG
    VerifyStiffnessSum(stiffness);
#endif
}
 
void TPZStructMatrixLCC::VerifyStiffnessSum(TPZMatrix<STATE> & stiffness){
    REAL totalSum = 0;
    for(int irow =0; irow < stiffness.Rows();irow++){
        for (int icol = 0; icol < stiffness.Cols(); icol++)
            totalSum += abs(stiffness.GetVal( irow, icol));
    }
    std::cout << "totalSum stisffnes =" << totalSum <<std::endl;
}
    
void TPZStructMatrixLCC::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
    
}

TPZMatrix<STATE> * TPZStructMatrixLCC::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *stiff = Create();
    
    //int64_t neq = stiff->Rows();
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

int TPZStructMatrixLCC::GetNumberColors(){
    int64_t nelem = fMesh->NElements();
        if (fElVecColor.size() != nelem) DebugStop();
        int ncolor = -1;
        for (int iel=0; iel<nelem; iel++){
            if (ncolor < fElVecColor[iel])
                ncolor = fElVecColor[iel];
        }
        ncolor++;
    return ncolor;
}

void TPZStructMatrixLCC::MultiThread_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
    if (fShouldColor){
        if (fUsingTBB){
            AssemblingUsingTBBandColoring(stiffness,rhs);
            
        }else{
            AssemblingUsingOMPandColoring(stiffness,rhs);
        }
        
    }else{
        if (fUsingTBB){
            AssemblingUsingTBBbutNotColoring(stiffness,rhs);
            
        }else{
            AssemblingUsingOMPbutNotColoring(stiffness,rhs);
        }
    }
    
    #ifdef PZDEBUG
        VerifyStiffnessSum(stiffness);
    #endif
    }
    
void TPZStructMatrixLCC::AssemblingUsingTBBandColoring(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs ){
#ifndef USING_TBB
    DebugStop();
#endif

    int64_t nelem = fMesh->NElements();
    const int nthread = this->fNumThreads;
    
    TPZStructMatrixLCC::OrderElements();
    int ncolor = GetNumberColors();
    
    for (int icol=0; icol<ncolor; icol++){
            
        tbb::task_scheduler_init init(nthread); //dont work in computer of LABMEC
        tbb::parallel_for( tbb::blocked_range<int64_t>(0,nelem),
                          [&](tbb::blocked_range<int64_t> r){
        for (int64_t iel = r.begin(); iel < r.end(); iel++)
            {
                if (icol != fElVecColor[iel]) continue;
                TPZCompEl *el = fMesh->Element(iel);
                if (!el) continue;
                
                ComputingCalcstiffAndAssembling(stiffness,rhs,el);

            }
        });
    }
    }


void TPZStructMatrixLCC::AssemblingUsingOMPandColoring(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs ){
#ifdef USING_OMP
    
    int64_t nelem = fMesh->NElements();
    const int nthread = this->fNumThreads;
            
    TPZStructMatrixLCC::OrderElements();
    int ncolor = GetNumberColors();
    
    for (int icol=0; icol<ncolor; icol++){
            
        omp_set_num_threads(nthread);
        #pragma omp parallel for schedule(dynamic,1)
        for (int64_t iel = 0; iel < nelem; iel++){
        if (icol != fElVecColor[iel]) continue;
        TPZCompEl *el = fMesh->Element(iel);
        if (!el) continue;
        
            ComputingCalcstiffAndAssembling(stiffness,rhs,el);

        }
    }
#else
    DebugStop();
#endif
    }

void TPZStructMatrixLCC::AssemblingUsingTBBbutNotColoring(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs ){
#ifdef USING_TBB
    int64_t nelem = fMesh->NElements();
    const int nthread = this->fNumThreads;
    

        tbb::task_scheduler_init init(nthread); //dont work in computer of LABMEC
        tbb::parallel_for( tbb::blocked_range<int64_t>(0,nelem),
                          [&](tbb::blocked_range<int64_t> r){
        for (int64_t iel = r.begin(); iel < r.end(); iel++)
        {
            TPZCompEl *el = fMesh->Element(iel);
            if (!el) continue;
            
            ComputingCalcstiffAndAssembling(stiffness,rhs,el);

        }
        });
#elif
    DebugStop();
#endif
    
}

void TPZStructMatrixLCC::AssemblingUsingOMPbutNotColoring(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs ){
#ifdef USING_OMP

    int64_t nelem = fMesh->NElements();
    const int nthread = this->fNumThreads;
    
    omp_set_num_threads(nthread);
    #pragma omp parallel for schedule(dynamic,1)
    for (int64_t iel = 0; iel < nelem; iel++){
        {
            TPZCompEl *el = fMesh->Element(iel);
            if (!el) continue;
                        
            ComputingCalcstiffAndAssembling(stiffness,rhs,el);
        }
    }
#else
    DebugStop();
#endif
}




void TPZStructMatrixLCC::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads(numthreads);
    int itr;
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return;
        }
    }
    

    this->fCurrentIndex = 0;
    fElementCompleted = -1;
    fElementsComputed.Resize(fMesh->NElements());
    fElementsComputed.Fill(0);
    fSomeoneIsSleeping = 0;
    TPZManVector<ThreadData*> allthreaddata(numthreads);

    for(itr=0; itr<numthreads; itr++)
    {
        allthreaddata[itr] = new ThreadData(this, itr, rhs, fMaterialIds, guiInterface);
        ThreadData &threaddata = *allthreaddata[itr];
        threaddata.fElBlocked=&fElBlocked;
        threaddata.fElSequenceColor=&fElSequenceColor;
        threaddata.fElementCompleted = &fElementCompleted;
        threaddata.fComputedElements = &fElementsComputed;
        threaddata.fSomeoneIsSleeping = &fSomeoneIsSleeping;

        allthreads.emplace_back(ThreadData::ThreadWork, &threaddata);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreads[itr].join();
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        delete allthreaddata[itr];
    }


#ifdef PZ_LOG
    if(loggerCheck.isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
    
}

void TPZStructMatrixLCC::ComputingCalcstiffAndAssembling(TPZMatrix<STATE>& stiffness,TPZFMatrix<STATE> &rhs,TPZCompEl *el){
    
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
    el->CalcStiff(ek, ef);
    
    if(!el->HasDependency()) {
        ek.ComputeDestinationIndices();
        fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
        
        if (fShouldColor){
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFelNonAtomic(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        }else{
            stiffness.AddKelAtomic(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        }
        
    } else {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
        ek.ComputeDestinationIndices();
        fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
        if (fShouldColor){
            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFelNonAtomic(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        }else{
            stiffness.AddKelAtomic(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        }
    
    }
}

TPZStructMatrixLCC::ThreadData::ThreadData(TPZStructMatrixLCC *strmat, int seqnum, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(&mat), fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
    fCurrentIndex = &strmat->fCurrentIndex;
    
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

TPZStructMatrixLCC::ThreadData::ThreadData(TPZStructMatrixLCC *strmat, int seqnum,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat),fGuiInterface(guiInterface), fGlobMatrix(0), fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
    this->fCurrentIndex = &strmat->fCurrentIndex;
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

//#define DRY_RUN

void *TPZStructMatrixLCC::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
//    TPZStructMatrixLCC *strmat = data->fStruct;
    TPZVec<int64_t> &ComputedElements = *(data->fComputedElements);
    TPZVec<int64_t> &ElBlocked = *(data->fElBlocked);

    int &SomeoneIsSleeping = *(data->fSomeoneIsSleeping);
    
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    int64_t numelements = data->fElSequenceColor->size();
    int64_t index = data->fCurrentIndex->fetch_add(1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
        sout << "index = " << index << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef HUGEDEBUG
    data->fMutexAccessElement.lock();
    std::cout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
    std::cout << "index = " << index << std::endl;
    std::cout.flush();
    data->fMutexAccessElement.unlock();
#endif
    
    while (index < numelements)
    {
        
        int64_t iel = data->fElSequenceColor->operator[](index);
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Computing element " << index;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef PZ_LOG
        std::stringstream sout;
        sout << "Element " << index << " elapsed time ";
        TPZTimer timeforel(sout.str());
        timeforel.start();
#endif

        if (iel >= 0){
            TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
            el->CalcStiff(ek,ef);
#else
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
            intel->InitializeElementMatrix(ek,ef);
#endif

            if(guiInterface) if(guiInterface->AmIKilled()){
                break;
            }
            
#ifndef DRY_RUN
            if(!el->HasDependency()) {
                
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef PZ_LOG
                if(loggerel.isDebugEnabled())
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
            }
#endif
            
#ifdef PZ_LOG
            timeforel.stop();
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << timeforel.processName() <<  timeforel;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            int64_t localcompleted = *(data->fElementCompleted);
            bool localupdated = false;
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
                localupdated = true;
            }
            int64_t needscomputed = ElBlocked[index];
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                if (localupdated) {
                    sout << "Localcompleted updated without thread lock\n";
                }
                
                sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
#ifdef HUGEDEBUG
            data->fMutexAccessElement.lock();
            std::cout << "threadEK " << data->fThreadSeqNum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
            data->fMutexAccessElement.unlock();
#endif
            
            bool hadtowait = false;
            while (needscomputed > localcompleted) {
                // block the thread till the element needed has been assembled
              std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
                SomeoneIsSleeping = 1;
                hadtowait = true;
#ifdef HUGEDEBUG
                std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
                std::cout.flush();
#endif
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element " << index << " cannot be assembled - going to sleep";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                data->fConditionVar.wait(lock);
                
                localcompleted = *data->fElementCompleted;
                localupdated = false;
                while (ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "thread wakeup for element index " << index << std::endl;
                    if (localupdated) {
                        sout << "Localcompleted updated without thread lock\n";
                    }
                    
                    sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
            
#ifdef HUGEDEBUG
            if (hadtowait) {
              std::scoped_lock lock(data->fMutexAccessElement);
                std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " continuing\n";
            }
#endif
            
#ifndef DRY_RUN
            auto globMatrix =
                    dynamic_cast<TPZMatrix<STATE> *>(data->fGlobMatrix);
                auto globRhs =
                    dynamic_cast<TPZFMatrix<STATE> *>(data->fGlobRhs);
            if(globMatrix){
                // Assemble the matrix
                if(!ek.HasDependency())
                {
                    globMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                    globRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                else
                {
                    globMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                    globRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                
            }
#endif
            
            localupdated = false;
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
                localupdated = true;
            }
            if (localcompleted == index-1) {
                localcompleted++;
            }
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
                localupdated = true;
            }
            bool elementcompletedupdate = false;
            if (*data->fElementCompleted < localcompleted) {
//                std::cout << "Updating element completed " << localcompleted << std::endl;
                *data->fElementCompleted = localcompleted;
                elementcompletedupdate = true;
            }
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element " << index << " has been assembled ";
                if (localupdated) {
                    sout << "\nLocalcompleted updated without thread lock\n";
                }
                if (elementcompletedupdate) {
                    sout << "\nfElementCompleted updated to localcompleted\n";
                }
                sout << "local completed " << localcompleted;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            ComputedElements[index] = 1;
            if (SomeoneIsSleeping) {
              std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
#ifdef HUGEDEBUG
                std::cout << "threadEK " <<data->fThreadSeqNum <<  " Computed index " << index << " Waking up ElementsCompleted " << *data->fElementCompleted << std::endl;
                std::cout.flush();
#endif
                
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    LOGPZ_DEBUG(logger, "condition broadcast")
                }
#endif
                SomeoneIsSleeping = 0;
                data->fConditionVar.notify_all();
            }

        }
        else
        {
            std::cout << "the element in ElColorSequence is negative???\n";
            DebugStop();
        }
        index = data->fCurrentIndex->fetch_add(1);
        
    }
    // just make sure threads that were accidentally blocked get woken up
    std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
    bool localupdated = false;
    int64_t localcompleted = *(data->fElementCompleted);
    while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
        localcompleted++;
        localupdated = true;
    }
    if (localcompleted == index-1) {
        localcompleted++;
    }
    bool elementcompletedupdate = false;
    if (*data->fElementCompleted < localcompleted) {
        //                std::cout << "Updating element completed " << localcompleted << std::endl;
        *data->fElementCompleted = localcompleted;
        elementcompletedupdate = true;
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger, "Finishing up")
        if (localupdated) {
            LOGPZ_DEBUG(logger, "updated localcompleted")
        }
        if (elementcompletedupdate) {
            LOGPZ_DEBUG(logger, "updated fElementCompleted")
        }
        if (localupdated || elementcompletedupdate) {
            std::stringstream sout;
            sout << "localcompleted " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
        LOGPZ_DEBUG(logger, "finishing and condition broadcast")
    }
#endif
    data->fConditionVar.notify_all();
    SomeoneIsSleeping = 0;
    return 0;
}

void *TPZStructMatrixLCC::ThreadData::ThreadWorkResidual(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    //    TPZStructMatrixLCC *strmat = data->fStruct;
    TPZVec<int64_t> &ComputedElements = *(data->fComputedElements);
    TPZVec<int64_t> &ElBlocked = *(data->fElBlocked);
    
    int &SomeoneIsSleeping = *(data->fSomeoneIsSleeping);
    
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    int64_t numelements = data->fElSequenceColor->size();
    int64_t index = data->fCurrentIndex->fetch_add(1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
        sout << "index = " << index << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef HUGEDEBUG
    data->fMutexAccessElement.lock();
    std::cout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
    std::cout << "index = " << index << std::endl;
    std::cout.flush();
    data->fMutexAccessElement.unlock();
#endif
    
        while (index < numelements)
        {
            
            int64_t iel = data->fElSequenceColor->operator[](index);
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Computing element " << index;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
#ifdef PZ_LOG
            std::stringstream sout;
            sout << "Element " << index << " elapsed time ";
            TPZTimer timeforel(sout.str());
            timeforel.start();
#endif
            
            if (iel >= 0){
                TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
                el->CalcResidual(ef);
#else
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
                intel->InitializeElementMatrix(ef);
#endif
                if(guiInterface) if(guiInterface->AmIKilled()){
                    break;
                }
                
#ifndef DRY_RUN
                if(!el->HasDependency()) {
                    
                    ef.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                    
#ifdef PZ_LOG
                    if(loggerel.isDebugEnabled())
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
                }
#endif
                
#ifdef PZ_LOG
                timeforel.stop();
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << timeforel.processName() <<  timeforel;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                int64_t localcompleted = *(data->fElementCompleted);
                bool localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                int64_t needscomputed = ElBlocked[index];
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    if (localupdated) {
                        sout << "Localcompleted updated without thread lock\n";
                    }
                    
                    sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
#ifdef HUGEDEBUG
                data->fMutexAccessElement.lock();
                std::cout << "threadEK " << data->fThreadSeqNum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
                data->fMutexAccessElement.unlock();
#endif
                
                bool hadtowait = false;
#ifdef PZ_LOG
                if (logger.isInfoEnabled() && needscomputed > localcompleted)
                {
                    std::stringstream sout;
                    sout << "Element " << index << " cannot be assembled - going to sleep";
                    LOGPZ_INFO(logger, sout.str())
                }
#endif
                while (needscomputed > localcompleted) {
                    // block the thread till the element needed has been assembled
                    std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
                    SomeoneIsSleeping = 1;
                    hadtowait = true;
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
                    std::cout.flush();
#endif
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Element " << index << " cannot be assembled - going to sleep";
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    data->fConditionVar.wait(lock);
                    
                    localcompleted = *data->fElementCompleted;
                    localupdated = false;
                    while (ComputedElements[localcompleted+1] == 1) {
                        localcompleted++;
                        localupdated = true;
                    }
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "thread wakeup for element index " << index << std::endl;
                        if (localupdated) {
                            sout << "Localcompleted updated without thread lock\n";
                        }
                        
                        sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                

#ifdef PZ_LOG
                if (logger.isInfoEnabled() && hadtowait)
                {
                    std::stringstream sout;
                    sout << "thread wakeup for element index " << index << std::endl;
                    LOGPZ_INFO(logger, sout.str())
                }
#endif


#ifdef HUGEDEBUG
                if (hadtowait) {
                  std::scoped_lock lock (data->fMutexAccessElement);
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " continuing\n";
                }
#endif
                
#ifndef DRY_RUN
                auto globRhs =
                    dynamic_cast<TPZFMatrix<STATE> *>(data->fGlobRhs);
                if(globRhs){
                    // Assemble the matrix
                    if(!ef.HasDependency())
                    {
                        globRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
                    }
                    else
                    {
                        globRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
                    }
                    
                }
#endif
                
                localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                if (localcompleted == index-1) {
                    localcompleted++;
                }
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                bool elementcompletedupdate = false;
                if (*data->fElementCompleted < localcompleted) {
                    //                std::cout << "Updating element completed " << localcompleted << std::endl;
                    *data->fElementCompleted = localcompleted;
                    elementcompletedupdate = true;
                }
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element " << index << " has been assembled ";
                    if (localupdated) {
                        sout << "\nLocalcompleted updated without thread lock\n";
                    }
                    if (elementcompletedupdate) {
                        sout << "\nfElementCompleted updated to localcompleted\n";
                    }
                    sout << "local completed " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                ComputedElements[index] = 1;
                if (SomeoneIsSleeping) {
                  std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum <<  " Computed index " << index << " Waking up ElementsCompleted " << *data->fElementCompleted << std::endl;
                    std::cout.flush();
#endif
                    
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        LOGPZ_DEBUG(logger, "condition broadcast")
                    }
#endif
                    SomeoneIsSleeping = 0;
                    data->fConditionVar.notify_all();
                }
                
            }
            else
            {
                std::cout << "the element in ElColorSequence is negative???\n";
                DebugStop();
            }
            index = data->fCurrentIndex->fetch_add(1);
            
        }
    // just make sure threads that were accidentally blocked get woken up
    std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
    bool localupdated = false;
    int64_t localcompleted = *(data->fElementCompleted);
    while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
        localcompleted++;
        localupdated = true;
    }
    if (localcompleted == index-1) {
        localcompleted++;
    }
    bool elementcompletedupdate = false;
    if (*data->fElementCompleted < localcompleted) {
        //                std::cout << "Updating element completed " << localcompleted << std::endl;
        *data->fElementCompleted = localcompleted;
        elementcompletedupdate = true;
    }
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger, "Finishing up")
        if (localupdated) {
            LOGPZ_DEBUG(logger, "updated localcompleted")
        }
        if (elementcompletedupdate) {
            LOGPZ_DEBUG(logger, "updated fElementCompleted")
        }
        if (localupdated || elementcompletedupdate) {
            std::stringstream sout;
            sout << "localcompleted " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
        LOGPZ_DEBUG(logger, "finishing and condition broadcast")
    }
#endif
    data->fConditionVar.notify_all();
    SomeoneIsSleeping = 0;
    return 0;
}

//static bool CanAssemble(TPZStack<int64_t> &connectlist, TPZVec<int> &elContribute)
//{
//    for (int i = 0 ; i < connectlist.NElements() ; i++)
//    {
//        if (elContribute[connectlist[i]] >= 0){
//            return false;
//        }
//    }
//    return true;
//}

static void AssembleColor(int64_t el, TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute) {
    for (auto connect : connectlist) {
        elContribute[connect] = el;
    }
}

static int64_t WhoBlockedMe(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &elSequenceColorInv) {
    int64_t el = -1;
    for (auto connect : connectlist) {
        auto elBlocked = elContribute[connect];
        if (elBlocked == -1) continue;
        auto elBlockedIndex = elSequenceColorInv[elBlocked];
        if (el == -1) el = elBlockedIndex;
        if (elBlockedIndex > el) el = elBlockedIndex;
    }
    return el;
}

//static void RemoveEl(int el,TPZCompMesh *cmesh,TPZVec<int> &elContribute,int elSequence)
//{
//    TPZCompEl *cel = cmesh->ElementVec()[el];
//    if(!cel) DebugStop();
//    TPZStack<int64_t> connectlist;
//    cel->BuildConnectList(connectlist);
//    for (int i = 0 ; i < connectlist.NElements() ; i++)
//    {
//        int conindex = connectlist[i];
//        if (elContribute[conindex] != elSequence){
//            DebugStop();
//        }
//        elContribute[conindex] = -1;
//    }
//}

static int64_t MinPassIndex(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &passIndex) {
    int64_t minPassIndex = -1;
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

//Coloring computational elements 2022
void TPZStructMatrixLCC::OrderElements(){
    int64_t nconnects = fMesh->NConnects();
    int64_t nelement = fMesh->ElementVec().NElements();
    fElVecColor.Resize(nelement);
    fElVecColor.Fill(-1);
    TPZVec<int> used(nconnects,0);
    bool haswork = true;
    int color = 0;
    
    while (haswork){
        haswork = false;
        for (int64_t iel = 0; iel < nelement; iel++){
            TPZCompEl *compel = fMesh->Element(iel);
            if (!compel) continue;
            if (fElVecColor[iel] != -1) continue;
            
            TPZStack<int64_t> connectlist;
            compel->BuildConnectList(connectlist);
            
            bool canColor = true;
            for (auto ic : connectlist){
                if (used[ic]){
                    canColor = false;
                    haswork = true;
                }
            }
            if (canColor){
                for (auto ic : connectlist) used[ic] = 1;
                fElVecColor[iel] = color;
            }
        }
        if (haswork) {
            color++;
            used.Fill(0);
        }
    }
}
 


int TPZStructMatrixLCC::ClassId() const{
    return Hash("TPZStructMatrixLCC") ^ TPZStructMatrixBase::ClassId() << 1;
}


void TPZStructMatrixLCC::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf, context);
    
    buf.Read(fElBlocked);
    buf.Read(fElSequenceColor);
    buf.Read(&fElementCompleted);
    buf.Read(&fSomeoneIsSleeping);
    int64_t fCurrentIndexLong;
    buf.Read(&fCurrentIndexLong);
    fCurrentIndex = fCurrentIndexLong;
}

void TPZStructMatrixLCC::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);
    
    buf.Write(fElBlocked);
    buf.Write(fElSequenceColor);
    buf.Write(&fElementCompleted);
    buf.Write(&fSomeoneIsSleeping);
    int64_t fCurrentIndexLong;
    fCurrentIndexLong = fCurrentIndex;
    buf.Write(&fCurrentIndexLong);
}

template class TPZRestoreClass<TPZStructMatrixLCC>;
