/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixST methods.
 */

#include "pzstrmatrixst.h"

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
#include "run_stats_table.h"

using namespace std;

#include "pzlog.h"

#include "pz_pthread.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixST"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


TPZStructMatrixST::TPZStructMatrixST(TPZCompMesh *mesh) : fMesh(mesh), fEquationFilter(mesh->NEquations()) {
    fMesh = mesh;
    this->SetNumThreads(0);
}

TPZStructMatrixST::TPZStructMatrixST(TPZAutoPointer<TPZCompMesh> cmesh) : fCompMesh(cmesh), fEquationFilter(cmesh->NEquations()) {
    fMesh = cmesh.operator->();
    this->SetNumThreads(0);
}

TPZStructMatrixST::TPZStructMatrixST(const TPZStructMatrixST &copy) : fMesh(copy.fMesh), fEquationFilter(copy.fEquationFilter)
{
    if (copy.fCompMesh) {
        fCompMesh = copy.fCompMesh;
    }
    fMaterialIds = copy.fMaterialIds;
    fNumThreads = copy.fNumThreads;
}

TPZMatrix<STATE> *TPZStructMatrixST::Create() {
    cout << "TPZStructMatrixST::Create should never be called\n";
    return 0;
}

TPZStructMatrixST *TPZStructMatrixST::Clone() {
    cout << "TPZStructMatrixST::Clone should never be called\n";
    return 0;
}


RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixST::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        long neqcondense = fEquationFilter.NActiveEquations();
#ifdef DEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense,rhs.Cols(),0.);

        this->OnlyAssemble(&stiffness, &rhsloc, guiInterface);
        
        fEquationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        this->OnlyAssemble(&stiffness, &rhs, guiInterface);
    }
    ass_stiff.stop();
}

void TPZStructMatrixST::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
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
        
        OnlyAssemble(&rhsloc, guiInterface);
        
        fEquationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
       OnlyAssemble(&rhs, guiInterface);
    }
    ass_rhs.stop();
}


void TPZStructMatrixST::ExecuteAssemble(TPZMatrix<STATE> *fGlobMatrix, TPZFMatrix<STATE> *fGlobRhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    
    
    TPZElementMatrix ek(fMesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(fMesh,TPZElementMatrix::EF);

    
    for(int iel=0; iel<fMesh->NElements(); iel++)
    {
        
            TPZCompEl *el = fMesh->ElementVec()[iel];
            
            if(!el) continue;
            
            if (fGlobMatrix)
                el->CalcStiff(ek,ef);
            else
                el->CalcResidual(ef);
            
            if(!el->HasDependency()) {
                
                if (fGlobMatrix) {
                    ek.ComputeDestinationIndices();
                    this->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                } else {
                    ef.ComputeDestinationIndices();
                    this->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                }
                
            } else {
                // the element has dependent nodes
                if (fGlobMatrix) {
                    ek.ApplyConstraints();
                    ef.ApplyConstraints();
                    ek.ComputeDestinationIndices();
                    this->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                } else {
                    ef.ApplyConstraints();
                    ef.ComputeDestinationIndices();
                    this->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                }
                
            }
            
            
            if(fGlobMatrix) {
                // assemble the matrix
                if(!ek.HasDependency()) {
                    fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                    fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                } else {
                    fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                    fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
            } else {
                if(!ef.HasDependency()) {
                    fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
                } else {
                    fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
                }
            }


        
    }
    
}


void TPZStructMatrixST::OnlyAssemble(TPZMatrix<STATE> *stiffness, TPZFMatrix<STATE> *rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ) {
    // Checking if the interface still exists
    if(guiInterface && guiInterface->AmIKilled()){
        return;
    }
    ExecuteAssemble(stiffness, rhs, guiInterface);
}

void TPZStructMatrixST::OnlyAssemble(TPZFMatrix<STATE> *rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    // Checking if the interface still exists
    if(guiInterface && guiInterface->AmIKilled()){
        return;
    }
    ExecuteAssemble(0, rhs, guiInterface);
}

/** Filter out the equations which are out of the range */
void TPZStructMatrixST::FilterEquations(TPZVec<long> &origindex, TPZVec<long> &destindex) const {
    fEquationFilter.Filter(origindex, destindex);
}

TPZMatrix<STATE> * TPZStructMatrixST::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    TPZMatrix<STATE> *stiff = Create();
    long cols = MAX(1, rhs.Cols());
    rhs.Redim(fEquationFilter.NEqExpand(),cols);
    Assemble(*stiff,rhs,guiInterface);
#ifdef LOG4CXX2
    if(loggerel->isDebugEnabled()) {
        std::stringstream sout;
        stiff->Print("Stiffness matrix",sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
#endif
    return stiff;
    
}

/** Set the set of material ids which will be considered when assembling the system */
void TPZStructMatrixST::SetMaterialIds(const std::set<int> &materialids) {
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


