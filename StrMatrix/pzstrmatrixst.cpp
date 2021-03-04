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


#include "pzcheckconsistency.h"
#include "TPZMaterial.h"
#include "run_stats_table.h"

using namespace std;

#include "pzlog.h"


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


TPZStructMatrixST::TPZStructMatrixST(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {
    
}

TPZStructMatrixST::TPZStructMatrixST(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
    
}

TPZStructMatrixST::TPZStructMatrixST(const TPZStructMatrixST &copy) : TPZStructMatrixBase(copy) {
    
}

TPZMatrix<STATE> *TPZStructMatrixST::Create() {
    cout << "TPZStructMatrixST::Create should never be called\n";
    return 0;
}

TPZStructMatrixST *TPZStructMatrixST::Clone() {
    cout << "TPZStructMatrixST::Clone should never be called\n";
    return 0;
}


static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixST::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
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
        int64_t neqcondense = fEquationFilter.NActiveEquations();
        int64_t neqexpand = fEquationFilter.NEqExpand();
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

#pragma omp parallel for private(ek, ef) schedule(dynamic, 1)
    for(int iel=0; iel<fMesh->NElements(); iel++)
	{
	
		TPZCompEl *el = fMesh->ElementVec()[iel];

		if(el) {

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

TPZMatrix<STATE> * TPZStructMatrixST::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    TPZMatrix<STATE> *stiff = Create();
    int64_t cols = MAX(1, rhs.Cols());
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

int TPZStructMatrixST::ClassId() const{
    return Hash("TPZStructMatrixST") ^ TPZStructMatrixBase::ClassId() << 1;
}

void TPZStructMatrixST::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf,context);
}

void TPZStructMatrixST::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);
}

template class TPZRestoreClass<TPZStructMatrixST>;