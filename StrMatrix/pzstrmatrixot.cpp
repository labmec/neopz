/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixOT methods.
 */

#include "pzstrmatrixot.h"
#include "pzsubcmesh.h"
#include "TPZElementMatrixT.h"
#include "TPZStructMatrix.h"

#include "TPZTimer.h"
#include "pzcheckconsistency.h"
#include "TPZMaterial.h"
#include "run_stats_table.h"

// #include <syncstream>
using namespace std;


#include <thread>
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixOT");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph_ot("-ass_graph_ot", "Run statistics table for the graph creation and coloring TPZStructMatrixOT.");


template<class TVar>
TPZStructMatrixOT<TVar>::TPZStructMatrixOT(): TPZStrMatParInterface(){}

template<class TVar>
TPZStructMatrixOT<TVar>::TPZStructMatrixOT(const TPZStructMatrixOT &copy) :
    TPZStrMatParInterface(copy),
    fElBlocked(copy.fElBlocked),
    fElSequenceColor(copy.fElSequenceColor),
    fElementsComputed(copy.fElementsComputed),
    fElementCompleted(copy.fElementCompleted)
    //fCurrentIndex //no need
{
}
//! Move constructor
template<class TVar>
TPZStructMatrixOT<TVar>::TPZStructMatrixOT(TPZStructMatrixOT &&copy) :
    TPZStrMatParInterface(copy),
    fElBlocked(copy.fElBlocked),
    fElSequenceColor(copy.fElSequenceColor),
    fElementsComputed(copy.fElementsComputed),
    fElementCompleted(copy.fElementCompleted)
    //fCurrentIndex //no need
{
}
//! Copy assignment operator
template<class TVar>
TPZStructMatrixOT<TVar>& TPZStructMatrixOT<TVar>::operator=(const TPZStructMatrixOT &copy){
    TPZStrMatParInterface::operator=(copy);
    fElBlocked = copy.fElBlocked;
    fElSequenceColor = copy.fElSequenceColor;
    fElementsComputed = copy.fElementsComputed;
    fElementCompleted = copy.fElementCompleted;
    //fCurrentIndex //no need
    return *this;
}
template<class TVar>
TPZStructMatrixOT<TVar>& TPZStructMatrixOT<TVar>::operator=(TPZStructMatrixOT &&copy){
    TPZStrMatParInterface::operator=(copy);
    fElBlocked = copy.fElBlocked;
    fElSequenceColor = copy.fElSequenceColor;
    fElementsComputed = copy.fElementsComputed;
    fElementCompleted = copy.fElementCompleted;
    //fCurrentIndex //no need
    return *this;
}


static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");


template<class TVar>
void TPZStructMatrixOT<TVar>::Assemble(TPZBaseMatrix & stiffness, TPZBaseMatrix & rhs){
    const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
    ass_stiff.start();
    if (equationFilter.IsActive()) {
        int64_t neqcondense = equationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<TVar> rhsloc;
        if(ComputeRhs()) rhsloc.Redim(neqcondense, rhs.Cols());
        if(this->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhsloc);
        }
        else{
            this->Serial_Assemble(stiffness,rhsloc);
        }
        
        if(ComputeRhs())equationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        if(this->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhs);
        }
        else{
            this->Serial_Assemble(stiffness,rhs);
        }
    }
    ass_stiff.stop();
}

template<class TVar>
void TPZStructMatrixOT<TVar>::Assemble(TPZBaseMatrix & rhs_base){
    const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
    if(!dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" Incompatible Types. Aborting...\n";
        DebugStop();
    }
    auto &rhs = dynamic_cast<TPZFMatrix<TVar> &>(rhs_base);
    ass_rhs.start();
    if(equationFilter.IsActive())
    {
        //TODONORM
        int64_t neqcondense = equationFilter.NActiveEquations();
        int64_t neqexpand = equationFilter.NEqExpand();
        if(rhs.Rows() != neqexpand || Norm(rhs) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<TVar> rhsloc(neqcondense,1,0.);
        if(this->fNumThreads)
        {
#ifdef HUGEDEBUG
            TPZFMatrix<TVar> rhsserial(rhsloc);
            this->Serial_Assemble(rhsserial);
#endif
            this->MultiThread_Assemble(rhsloc);
#ifdef HUGEDEBUG
            rhsserial -= rhsloc;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
        }
        else
        {
            this->Serial_Assemble(rhsloc);
        }
        equationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        if(this->fNumThreads){
#ifdef HUGEDEBUG
            TPZFMatrix<TVar> rhsserial(rhs);
            this->Serial_Assemble(rhsserial);
#endif
            this->MultiThread_Assemble(rhs);
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
            this->Serial_Assemble(rhs);
        }
    }
    ass_rhs.stop();
}

template<class TVar>
void TPZStructMatrixOT<TVar>::InitCreateAssemble()
{
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    stat_ass_graph_ot.start();
    TPZManVector<int64_t> ElementOrder;
    TPZStructMatrixOT<TVar>::OrderElement(myself->Mesh(), ElementOrder);
    TPZVec<int64_t> elcolors;
    TPZStructMatrixOT<TVar>::ElementColoring(myself->Mesh(), ElementOrder, fElSequenceColor, fElBlocked, elcolors);
    stat_ass_graph_ot.stop();
}

template<class TVar>
void TPZStructMatrixOT<TVar>::Serial_Assemble(TPZBaseMatrix & stiff_base, TPZBaseMatrix & rhs_base){
    if(!dynamic_cast<TPZMatrix<TVar>*>(&stiff_base)||
       !dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<": incompatible types. Aborting...\n";
        DebugStop();
    }
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();
    const auto &equationFilter = myself->EquationFilter();
    auto &materialIds = myself->MaterialIds();
    
    auto &stiffness = dynamic_cast<TPZMatrix<TVar>&>(stiff_base);
    auto &rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
    if(!cmesh){
        LOGPZ_ERROR(logger,"Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef PZ_LOG
    if(dynamic_cast<TPZSubCompMesh * >(cmesh))
    {
        std::stringstream sout;
        sout << "AllEig = {};";
        LOGPZ_DEBUG(loggerelmat,sout.str())
        
    }
#endif
#ifdef PZDEBUG
    if (ComputeRhs() && rhs.Rows() != equationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif
    
    int64_t iel;
    int64_t nelem = cmesh->NElements();
    TPZElementMatrixT<TVar> ek(cmesh, TPZElementMatrix::EK),ef(cmesh, TPZElementMatrix::EF);
#ifdef PZ_LOG
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    
    int64_t count = 0;
    for(iel=0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        int matidsize = materialIds.size();
        if(matidsize){
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat)
            {
                if (!submesh) {
                    continue;
                }
                else if(submesh->NeedsComputing(materialIds) == false) continue;
            }
            else
            {
                int matid = mat->Id();
                if (myself->ShouldCompute(matid) == false) continue;
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
        
        
#ifdef PZ_LOG
        if(dynamic_cast<TPZSubCompMesh * >(cmesh))
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
            equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //			TPZSFMatrix<TVar> test(stiffness);
            //			TPZFMatrix<TVar> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //			stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            if(ComputeRhs())rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //			rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //			test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			test -= stiffness;
            //			test.Print("matriz de rigidez diference",std::cout);
            //			test2.Print("matriz de rigidez interface",std::cout);
            
#ifdef PZ_LOG
            if(loggerel.isDebugEnabled())
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
                if(ComputeRhs())ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            if(ComputeRhs())rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef PZ_LOG
            if(loggerel.isDebugEnabled() && ! dynamic_cast<TPZSubCompMesh *>(cmesh))
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
    
#ifdef PZ_LOG
    if(loggerCheck.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        stiffness.Print("Matriz de Rigidez: ",sout,EMathematicaInput);
        if(ComputeRhs()) rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
    
#endif
    
}
template<class TVar>
void TPZStructMatrixOT<TVar>::Serial_Assemble(TPZBaseMatrix & rhs_base){
    if(!dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<": incompatible types. Aborting...\n";
        DebugStop();
    }

    auto &rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    const auto &equationFilter = myself->EquationFilter();
    auto *cmesh = myself->Mesh();
    int64_t iel;
//    int64_t nelem = cmesh->NElements();
    
    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");
    
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    
    stat_ass_graph_ot.start();
    TPZManVector<int64_t> ElementOrder;
    TPZStructMatrixOT<TVar>::OrderElement(cmesh, ElementOrder);
    TPZVec<int64_t> elcolors;
    TPZStructMatrixOT<TVar>::ElementColoring(cmesh, ElementOrder, fElSequenceColor, fElBlocked, elcolors);
    stat_ass_graph_ot.stop();

    int64_t elseqsize = fElSequenceColor.size();
    for (int64_t index = 0; index < elseqsize; index++) {
        iel = fElSequenceColor[index];
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        
        TPZMaterial * mat = el->Material();
        if (!mat) continue;
        int matid = mat->Id();
        if (myself->ShouldCompute(matid) == false) continue;
        
        TPZElementMatrixT<TVar> ef(cmesh, TPZElementMatrix::EF);
        
        calcresidual.start();
        
        el->CalcResidual(ef);
        
        calcresidual.stop();
        
        assemble.start();
        
        if(!el->HasDependency()) {
            ef.ComputeDestinationIndices();
            equationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            equationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
        }
        
        assemble.stop();
        
    }//fim for iel
#ifdef PZ_LOG
    {
        if(logger.isDebugEnabled())
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


template<class TVar>
void TPZStructMatrixOT<TVar>::MultiThread_Assemble(TPZBaseMatrix & mat_base, TPZBaseMatrix & rhs_base)
{
    if(!dynamic_cast<TPZMatrix<TVar>*>(&mat_base)||
       !dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<": incompatible types. Aborting...\n";
        DebugStop();
    }

    auto &mat = dynamic_cast<TPZMatrix<TVar>&>(mat_base);
    auto &rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);

    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    const int numthreads = this->fNumThreads;
    // std::cout << "Assemble numthreads = " << numthreads << std::endl;
    std::vector<std::thread> allthreads;
    int itr;

    this->fCurrentIndex = 0;
    
    fElementCompleted = -1;
    fElementsComputed.Resize(myself->Mesh()->NElements()+1);
    fElementsComputed.Fill(0);
#ifdef HUGEDEBUG
    {
        for (int64_t i=1; i<fElBlocked.size(); i++) {
            std::cout << "i = " << i << " fElBlocked[i] " << fElBlocked[i] << std::endl;
        }
    }
#endif
    TPZManVector<ThreadData *> allthreadsD(numthreads,0); 

    std::mutex accessElement;
    std::condition_variable conditionVar;
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreadsD[itr] = new ThreadData(myself, itr, mat, rhs, myself->MaterialIds(),
                          &fCurrentIndex, ComputeRhs());
        allthreadsD[itr]->fElBlocked= &fElBlocked;
        allthreadsD[itr]->fElSequenceColor=&fElSequenceColor;
        allthreadsD[itr]->fElementCompleted = &fElementCompleted;
        allthreadsD[itr]->fComputedElements = &fElementsComputed;
        allthreadsD[itr]->fMutexAccessElement = &accessElement;
        allthreadsD[itr]->fConditionVar = &conditionVar;

        allthreadsD[itr]->fThreadSeqNum = itr;
        allthreads.push_back(std::thread(ThreadData::ThreadWork, allthreadsD[itr]));
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreads[itr].join();
    }

#ifdef PZ_LOG
    if(loggerCheck.isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        mat.Print("Matrix: ",sout,EMathematicaInput);
        if(ComputeRhs())rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
}

template<class TVar>
void TPZStructMatrixOT<TVar>::MultiThread_Assemble(TPZBaseMatrix & rhs)
{
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads(numthreads);
    int itr;    

    this->fCurrentIndex = 0;
    fElementCompleted = -1;
    fElementsComputed.Resize(myself->Mesh()->NElements());
    fElementsComputed.Fill(0);
    TPZManVector<ThreadData*> allthreaddata(numthreads);

    for(itr=0; itr<numthreads; itr++)
    {
        allthreaddata[itr] =
            new ThreadData(myself, itr, rhs,myself->MaterialIds(),&fCurrentIndex);
        ThreadData &threaddata = *allthreaddata[itr];
        threaddata.fElBlocked=&fElBlocked;
        threaddata.fElSequenceColor=&fElSequenceColor;
        threaddata.fElementCompleted = &fElementCompleted;
        threaddata.fComputedElements = &fElementsComputed;

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



template<class TVar>
TPZStructMatrixOT<TVar>::ThreadData::ThreadData(
    TPZStructMatrix *strmat, int seqnum, TPZBaseMatrix &mat,
    TPZBaseMatrix &rhs,const std::set<int> &MaterialIds,
    std::atomic<int64_t> *currentIndex, bool computeRhs)
: fStruct(strmat), fGlobMatrix(&mat),
  fGlobRhs(&rhs), fThreadSeqNum(seqnum), fComputeRhs(computeRhs)
{
    fCurrentIndex = currentIndex;
    
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

template<class TVar>
TPZStructMatrixOT<TVar>::ThreadData::ThreadData(
    TPZStructMatrix *strmat, int seqnum, TPZBaseMatrix &rhs,
    const std::set<int> &MaterialIds,
    std::atomic<int64_t> *currentIndex)
: fStruct(strmat),fGlobMatrix(0),
  fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
    this->fCurrentIndex = currentIndex;
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

#ifdef USING_MKL
#include <mkl.h>
#endif
//#define HUGEDEBUG

template<class TVar>
void *TPZStructMatrixOT<TVar>::ThreadData::ThreadWork(void *datavoid)
{
#ifdef USING_MKL
    //we keep the current value of n threads
    auto mkl_threads = mkl_set_num_threads_local(1);
#endif
    ThreadData *data = (ThreadData *) datavoid;
    int seqnum = data->fThreadSeqNum;
    TPZVec<int64_t> &ComputedElements = *(data->fComputedElements);
    TPZVec<int64_t> &ElBlocked = *(data->fElBlocked);

    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZElementMatrixT<TVar> ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrixT<TVar> ef(cmesh,TPZElementMatrix::EF);
    int64_t numelements = data->fElSequenceColor->size();
    int64_t index = data->fCurrentIndex->fetch_add(1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Thread number " << seqnum << " total elements " << numelements;
        sout << " locked index = " << index;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef HUGEDEBUG
    {
        std::stringstream sout;
        sout << "Thread " << seqnum << " total elements " << numelements << "index = " << index << std::endl;
        std::cout << sout.str();
        // std::osyncstream(std::cout) << sout.str();
    }
#endif

    while (index < numelements)
    {
        // {
        //     sout << "Thread " << seqnum << " computing index " << index << std::endl;
        //     sout.flush();
        // }
        int64_t iel = data->fElSequenceColor->operator[](index);
        if(iel < 0) DebugStop();
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Computing element " << index;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
        el->CalcStiff(ek,ef);
#else
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
        intel->InitializeElementMatrix(ek,ef);
#endif
        
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
        

        int64_t localcompleted = *(data->fElementCompleted);
        bool localupdated = false;
        while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
            localcompleted++;
            localupdated = true;
        }
        if(localcompleted > *(data->fElementCompleted)) {
            *(data->fElementCompleted) = localcompleted;
        }
        const int64_t index_required = ElBlocked[index];
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "thread " << seqnum << " localupdated " << localupdated;
            sout << " Element " << index << " is computed and can assemble if required " << index_required << " is smaller than localcompleted " << localcompleted << std::endl;
            std::cout << sout.str();
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
#ifdef HUGEDEBUG
    {
        std::stringstream sout;
        sout << "threadEK " << seqnum << " index " << index << " localcompleted " << localcompleted << " index required " << index_required << std::endl;
        std::cout << sout.str();
        // std::osyncstream(std::cout) << "threadEK " << seqnum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
    }
#endif
        
        if (index_required > localcompleted) {
#ifdef HUGEDEBUG
            {
                std::stringstream sout; sout << "Thread " << seqnum << " trying to acquire lock\n";
                cout << sout.str();
            }
#endif
            // block the thread till the element needed has been assembled
            std::unique_lock<std::mutex> lock(*data->fMutexAccessElement);
#ifdef HUGEDEBUG
            std::cout << "threadEK " << seqnum << " Index " << index << " going to sleep waiting for " << index_required << std::endl;
            // std::osyncstream(std::cout) << "threadEK " << seqnum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
#endif
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element " << index << " cannot be assembled - going to sleep";
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
#ifdef HUGEDEBUG
            {
                std::stringstream sout; sout << "Thread " << seqnum << " acquired lock\n";
                std::cout << sout.str();
            }
#endif

            data->fConditionVar->wait(lock,[index_required,&data](){return index_required <=*(data->fElementCompleted);});
            
        }
#ifdef HUGEDEBUG
        {
            std::stringstream sout; sout << "Thread " << seqnum << " woke up\n";
            std::cout << sout.str();
        }
#endif

        
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "thread " << seqnum << " wakeup element " << index << std::endl;
            sout << "Element " << index << " is computed and can assemble if required " << index_required << " is smaller than element completed " << (*data->fElementCompleted) << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
        
#ifndef DRY_RUN
        auto globMatrix =
                dynamic_cast<TPZMatrix<TVar> *>(data->fGlobMatrix);
            auto globRhs =
                dynamic_cast<TPZFMatrix<TVar> *>(data->fGlobRhs);
        if(globMatrix){
            // Assemble the matrix
            if(!ek.HasDependency())
            {
                globMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                if(data->fComputeRhs) globRhs->AddFelNonAtomic(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            }
            else
            {
                globMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                if(data->fComputeRhs)globRhs->AddFelNonAtomic(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            }
            
        }
#endif
        ComputedElements[index] = 1;

        {
            int64_t localcompleted = *(data->fElementCompleted);
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
            }
            *(data->fElementCompleted) = localcompleted;
        }
        data->fConditionVar->notify_all();
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Element " << index << " has been assembled ";
            sout << "local completed " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef HUGEDEBUG
        {
            // std::osyncstream(std::cout) << "Thread " << seqnum << " Element " << index << " has been assembled\n";
            std::stringstream sout; sout << "Thread " << seqnum << " Element " << index << " has been assembled fElementCompleted " << *(data->fElementCompleted) << "\n";
            std::cout << sout.str();
            // sout.flush();
        }
#endif

        index = data->fCurrentIndex->fetch_add(1);
#ifdef HUGEDEBUG
        {
            std::stringstream sout;
            sout << "Thread " << seqnum << " next index = " << index << " numelements " << numelements << std::endl;
            std::cout << sout.str();
        }
#endif        
    } // while (index < numelements)
    // just make sure threads that were accidentally blocked get woken up
#ifdef HUGEDEBUG
    {
        std::stringstream sout; sout << "Thread " << seqnum << " trying to acquire lock " << std::endl;
        std::cout << sout.str();
    }
#endif
    {
        std::scoped_lock<std::mutex> lock(*data->fMutexAccessElement);
    
#ifdef HUGEDEBUG
        {
            std::stringstream sout; sout << "Thread " << seqnum << " acquired lock " << std::endl;
            std::cout << sout.str();
        }
#endif
        int64_t localcompleted = *(data->fElementCompleted);
        while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
            localcompleted++;
        }
        if (*data->fElementCompleted < localcompleted) {
#ifdef HUGEDEBUG
            std::cout << "Updating element completed " << localcompleted << std::endl;
#endif
            *data->fElementCompleted = localcompleted;
        }

#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            LOGPZ_DEBUG(logger, "Finishing up")
            LOGPZ_DEBUG(logger, "finishing and condition broadcast")
        }
#endif
#ifdef HUGEDEBUG
        std::cout << "Thread " << seqnum << " Bailing out index = " << index << " numelements " << numelements << std::endl;
#endif
        // std::osyncstream(std::cout) << "Thread " << seqnum << " Bailing out index = " << index << " numelements " << numelements << std::endl;
        data->fConditionVar->notify_all();
#ifdef HUGEDEBUG
        std::cout << "Thread " << seqnum << " releasing lock\n";
#endif
    }
    
#ifdef USING_MKL
    //we restore the previous value of n threads
    mkl_set_num_threads_local(mkl_threads);
#endif
    return 0;
}


template<class TVar>
TPZStructMatrixOT<TVar>::ThreadData::~ThreadData() {
}

template<class TVar>
bool TPZStructMatrixOT<TVar>::ThreadData::ShouldCompute(int matid) const{
    return fStruct->ShouldCompute(matid);
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

/**
 * elSequence (input) element sequence according to the connect sequence numbers
 * elSequenceColor (output) the colored element sequence
 * elBlocked the element index which needs to have been computed before assembling the element
 * elColors (output) number of elements in each color
 */
template<class TVar>
void TPZStructMatrixOT<TVar>::ElementColoring(TPZCompMesh *cmesh, TPZVec<int64_t> &elSequence, TPZVec<int64_t> &elSequenceColor,
        TPZVec<int64_t> &elBlocked, TPZVec<int64_t> &NumelColors) {

    const int64_t nnodes = cmesh->NConnects();
    const int64_t nel = cmesh->NElements();

    if (nel == 0) return;

    NumelColors.Resize(nel, -1);

    // elContribute contains the element index which last contributed to the node
    // passIndex contains the color of the element
    TPZManVector<int64_t> elContribute(nnodes, -1), passIndex(nel, -1), elSequenceColorInv(nel, -1);
    int64_t elsequencesize = elSequence.size();
    elSequenceColor.Resize(elsequencesize);
    elSequenceColor.Fill(-1);
    elBlocked.Resize(elsequencesize);
    elBlocked.Fill(-1);
    int64_t nelProcessed = 0;
    int64_t currentPassIndex = 0;
    while (nelProcessed < elSequence.NElements()) {
        for (auto elindex : elSequence) {
            // if this element hasn't been computed in a previous pass
            if (elSequenceColorInv[elindex] == -1) {
                TPZCompEl *cel = cmesh->Element(elindex);
                if (!cel) continue;
                TPZStack<int64_t> connectlist;
                cel->BuildConnectList(connectlist);
                // compute the lowest color (pass index) of the elements that have contributed to this set of nodes
                int64_t minPass = MinPassIndex(connectlist, elContribute, passIndex);
                // no element has ever seen any of these nodes
                if (minPass == -1) {
                    passIndex[elindex] = currentPassIndex;
                    // the element index is put into elContribute (from there it is possible to know the color as well)
                    AssembleColor(elindex, connectlist, elContribute);
                    // initialize the data structures
                    elSequenceColor[nelProcessed] = elindex;
                    elSequenceColorInv[elindex] = nelProcessed;
                    nelProcessed++;
                } // this element cannot be computed as it is connected to another element of the current color
                else if (minPass == currentPassIndex) {
                } // the elements connected to this node are from a previous color
                else if (minPass < currentPassIndex) {
                    // the element with largest index which contributes to the set of nodes
                    // el is given in the new sequence order
                    const int64_t el = WhoBlockedMe(connectlist, elContribute, elSequenceColorInv);
                    // elblocked means the future element the element el will block
                    if (elBlocked[nelProcessed] != -1) DebugStop();
                    elBlocked[nelProcessed] = el;
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
        NumelColors[currentPassIndex] = nelProcessed;
        currentPassIndex++;
    }

    NumelColors[currentPassIndex] = NumelColors[currentPassIndex - 1] + 1;
    NumelColors.Resize(currentPassIndex + 1);

// #ifdef PZDEBUG
//     std::ofstream toto("../ColorMeshDebug.txt");
//     toto << "elSequence\n" << elSequence << std::endl;
//     toto << "elSequenceColor\n" << elSequenceColor << std::endl;
//     toto << "elSequenceColorInv\n" << elSequenceColorInv << std::endl;
//     toto << "elBlocked\n" << elBlocked << std::endl;
//     toto << "elContribute\n" << elContribute << std::endl;
//     toto << "passIndex\n" << passIndex << std::endl;
//     toto << "NumelColors\n" << NumelColors << std::endl;
//     toto.close();
// #endif
}

template<class TVar>
void TPZStructMatrixOT<TVar>::OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder) {
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
    int64_t el;
    TPZCompEl *cel;

#ifdef NOORDER
    int64_t count = 0;
    std::cout << __PRETTY_FUNCTION__ << " no element order\n";
    ElementOrder.Resize(cmesh->NElements(), -1);
    for (el = 0; el < cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        ElementOrder[count] = el;
        count++;
    }
    ElementOrder.Resize(count);
    return;
#endif
    for (el = 0; el < cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
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
    TPZVec<int> elorderinv(cmesh->ElementVec().NElements(), -1);
    for (int64_t seq = 0; seq < nconnect; seq++) {
        int64_t ic = nodeorder[seq];
        if (ic == -1) continue;
        int64_t firstind = firstelconnect[ic];
        int64_t lastind = firstelconnect[ic + 1];
        for (int64_t ind = firstind; ind < lastind; ind++) {
            el = elconnect[ind];
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
template<class TVar>
int TPZStructMatrixOT<TVar>::ClassId() const{
    return Hash("TPZStructMatrixOT") ^
        TPZStrMatParInterface::ClassId() << 1 ^
        ClassIdOrHash<TVar>() << 2;
}

template<class TVar>
void TPZStructMatrixOT<TVar>::Read(TPZStream& buf, void* context) {
    TPZStrMatParInterface::Read(buf, context);
    
    buf.Read(fElBlocked);
    buf.Read(fElSequenceColor);
    buf.Read(&fElementCompleted);
    int64_t fCurrentIndexLong;
    buf.Read(&fCurrentIndexLong);
    fCurrentIndex = fCurrentIndexLong;
}

template<class TVar>
void TPZStructMatrixOT<TVar>::Write(TPZStream& buf, int withclassid) const {
    TPZStrMatParInterface::Write(buf, withclassid);
    
    buf.Write(fElBlocked);
    buf.Write(fElSequenceColor);
    buf.Write(&fElementCompleted);
    int64_t fCurrentIndexLong;
    fCurrentIndexLong = fCurrentIndex;
    buf.Write(&fCurrentIndexLong);
}

template class TPZStructMatrixOT<STATE>;
template class TPZRestoreClass<TPZStructMatrixOT<STATE>>;
template class TPZStructMatrixOT<CSTATE>;
template class TPZRestoreClass<TPZStructMatrixOT<CSTATE>>;
