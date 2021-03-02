/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixCS methods.
 */

#include "pzstrmatrixcs.h"

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
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixCS"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


TPZStructMatrixCS::TPZStructMatrixCS(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {
    
}

TPZStructMatrixCS::TPZStructMatrixCS(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
    
}

TPZStructMatrixCS::TPZStructMatrixCS(const TPZStructMatrixCS &copy) : TPZStructMatrixBase(copy) {
    
}

TPZMatrix<STATE> *TPZStructMatrixCS::Create() {
    cout << "TPZStructMatrixCS::Create should never be called\n";
    return 0;
}

TPZStructMatrixCS *TPZStructMatrixCS::Clone() {
    cout << "TPZStructMatrixCS::Clone should never be called\n";
    return 0;
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixCS::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){

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

void TPZStructMatrixCS::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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



void TPZStructMatrixCS::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ){
    
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

void TPZStructMatrixCS::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
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

TPZMatrix<STATE> * TPZStructMatrixCS::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
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


void TPZStructMatrixCS::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    ThreadData threaddata(this,mat,rhs,fMaterialIds,guiInterface);
#ifdef USING_TBB
    AssembleTask tasks(&threaddata);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, fMesh->NElements()), tasks);
#else
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
#endif
    
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


void TPZStructMatrixCS::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    ThreadData threaddata(this,rhs,fMaterialIds,guiInterface);
    
#ifdef USING_TBB
    AssembleTask tasks(&threaddata);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, fMesh->NElements()), tasks);
#else
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
#endif
}



TPZStructMatrixCS::ThreadData::ThreadData(TPZStructMatrixCS *strmat, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0)
{

}

TPZStructMatrixCS::ThreadData::ThreadData(TPZStructMatrixCS *strmat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat),fGuiInterface(guiInterface), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0)
{
  
}


void *TPZStructMatrixCS::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    // compute the next element (this method is threadsafe)
    int64_t iel = data->NextElement();
    TPZCompMesh *cmesh = data->fStruct->fMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    int64_t nel = cmesh->NElements();
    while(iel < nel)
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

        TPZAutoPointer<TPZElementMatrix> ek;
        TPZAutoPointer<TPZElementMatrix> ef = new TPZElementMatrix(cmesh,TPZElementMatrix::EF);
        if (data->fGlobMatrix) {
            ek = new TPZElementMatrix(cmesh,TPZElementMatrix::EK);
        }
        else
        {
            ek = ef;
        }
        
        TPZCompEl *el = cmesh->ElementVec()[iel];
        TPZElementMatrix *ekp = ek.operator->();
        TPZElementMatrix *efp = ef.operator->();
        TPZElementMatrix &ekr = *ekp;
        TPZElementMatrix &efr = *efp;
        
        if (data->fGlobMatrix) {
            el->CalcStiff(ekr,efr);
        }
        else
        {
            el->CalcResidual(efr);
        }
        
        if(guiInterface) if(guiInterface->AmIKilled()){
            break;
        }
        
        if(!el->HasDependency()) {
            ek->ComputeDestinationIndices();
            
            if(data->fStruct->HasRange())
            {
                data->fStruct->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
            }
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element index " << iel << std::endl;
                ek->fMat.Print("Element stiffness matrix",sout);
                ef->fMat.Print("Element right hand side", sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            if (data->fGlobMatrix) {
                ek->ApplyConstraints();
            }
            ef->ApplyConstraints();
            ek->ComputeDestinationIndices();
            if(data->fStruct->HasRange())
            {
                data->fStruct->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
            }
#ifdef LOG4CXX
            if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
            {
                std::stringstream sout;
                el->Reference()->Print(sout);
                el->Print(sout);
                ek->Print(sout);
                //			ef->Print(sout);
                LOGPZ_DEBUG(loggerel2,sout.str())
            }
#endif
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element index " << iel << std::endl;
                ek->fConstrMat.Print("Element stiffness matrix",sout);
                ef->fConstrMat.Print("Element right hand side", sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        }
        
        // Assemble the matrix
        
        if(!ek->HasDependency())
        {
            if (data->fGlobMatrix) {
              std::scoped_lock lock(data->fMutexAccessElementK);
                data->fGlobMatrix->AddKel(ek->fMat,ek->fSourceIndex,ek->fDestinationIndex);                
            }
            std::scoped_lock lock(data->fMutexAccessElementF);
            data->fGlobRhs->AddFel(ef->fMat,ek->fSourceIndex,ek->fDestinationIndex);
        }
        else
        {
            if (data->fGlobMatrix) {
                std::scoped_lock lock(data->fMutexAccessElementK);
                data->fGlobMatrix->AddKel(ek->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
            }
            std::scoped_lock lock(data->fMutexAccessElementF);
            data->fGlobRhs->AddFel(ef->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
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

        // compute the next element (this method is threadsafe)
        iel = data->NextElement();
    }
    
    return 0;
}

int64_t TPZStructMatrixCS::ThreadData::NextElement()
{
    TPZCompMesh *cmesh = fStruct->fMesh;
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    
    int64_t nel = elementvec.NElements();
    int64_t my_el;
    
    while (1) {
      {
        std::scoped_lock lock(fMutexAccessElement);
        my_el = fNextElement++;
      }
        
        
        if (my_el >= nel-1)
            break;
        
        TPZCompEl *el = elementvec[my_el];
        if(!el) continue;
        if(fStruct->fMaterialIds.size() == 0) break;
        TPZMaterial * mat = el->Material();
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
        if(!mat)
        {
            if(!submesh)
            {
                continue;
            }
            else if(submesh->NeedsComputing(fStruct->fMaterialIds) == false) continue;
        }
        else
        {
            int matid = mat->Id();
            if(this->ShouldCompute(matid) == false) continue;
        }
        break;
    }
    
    return my_el;
}

#ifdef USING_TBB
void TPZStructMatrixCS::AssembleTask::operator()(const tbb::blocked_range<size_t>& range) const
{
    TPZCompMesh *cmesh = data->fStruct->fMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    
    for(size_t iel=range.begin(); iel!=range.end(); ++iel )
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

        TPZAutoPointer<TPZElementMatrix> ek;
        TPZAutoPointer<TPZElementMatrix> ef = new TPZElementMatrix(cmesh,TPZElementMatrix::EF);
        if (data->fGlobMatrix) {
            ek = new TPZElementMatrix(cmesh,TPZElementMatrix::EK);
        }
        else
        {
            ek = ef;
        }
        
        TPZCompEl *el = cmesh->ElementVec()[iel];
        
        if (el) {
            
        TPZElementMatrix *ekp = ek.operator->();
        TPZElementMatrix *efp = ef.operator->();
        TPZElementMatrix &ekr = *ekp;
        TPZElementMatrix &efr = *efp;
        
        if (data->fGlobMatrix) {
            el->CalcStiff(ekr,efr);
        }
        else
        {
            el->CalcResidual(efr);
        }
        
        if(guiInterface) if(guiInterface->AmIKilled()){
            break;
        }
        
        if(!el->HasDependency()) {
            ek->ComputeDestinationIndices();
            
            if(data->fStruct->HasRange())
            {
                data->fStruct->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
            }
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element index " << iel << std::endl;
                ek->fMat.Print("Element stiffness matrix",sout);
                ef->fMat.Print("Element right hand side", sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            if (data->fGlobMatrix) {
                ek->ApplyConstraints();
            }
            ef->ApplyConstraints();
            ek->ComputeDestinationIndices();
            if(data->fStruct->HasRange())
            {
                data->fStruct->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
            }
#ifdef LOG4CXX
            if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
            {
                std::stringstream sout;
                el->Reference()->Print(sout);
                el->Print(sout);
                ek->Print(sout);
                //			ef->Print(sout);
                LOGPZ_DEBUG(loggerel2,sout.str())
            }
#endif
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element index " << iel << std::endl;
                ek->fConstrMat.Print("Element stiffness matrix",sout);
                ef->fConstrMat.Print("Element right hand side", sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
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
        // Assemble the matrix
        
        if(!ek->HasDependency())
        {
            if (data->fGlobMatrix) {
              data->fGlobMatrix->AddKel(ek->fMat,ek->fSourceIndex,ek->fDestinationIndex);
                                
            }
            data->fGlobRhs->AddFel(ef->fMat,ek->fSourceIndex,ek->fDestinationIndex);
        }
        else
        {
            if (data->fGlobMatrix) {
              data->fGlobMatrix->AddKel(ek->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
            }
            
            data->fGlobRhs->AddFel(ef->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
            
        }
        }
    }
}
#endif

int TPZStructMatrixCS::ClassId() const{
    return Hash("TPZStructMatrixCS") ^ TPZStructMatrixBase::ClassId() << 1;
}

void TPZStructMatrixCS::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf, context);
    
}

void TPZStructMatrixCS::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);

}

template class TPZRestoreClass<TPZStructMatrixCS>;
