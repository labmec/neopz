/**
 * @file
 * @brief Contains implementations of the TPZSubMeshFrontalAnalysis methods.
 */

#include "pzsmfrontalanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
#include "TPZMatrixSolver.h"

using namespace std;

// Construction/Destruction

TPZSubMeshFrontalAnalysis::TPZSubMeshFrontalAnalysis(TPZSubCompMesh *mesh) : TPZLinearAnalysis(mesh)                                                , fReferenceSolution(false)//TODOCOMPLEX:set matrix type (complex/real)
{
	fMesh = mesh;
	fFront = 0;
}

TPZSubMeshFrontalAnalysis::~TPZSubMeshFrontalAnalysis()
{
	
}

void TPZSubMeshFrontalAnalysis::Run(std::ostream &out){
	
	//fReducableStiff.Print("Reducable stiff before assembled");
    //TODOCOMPLEX
	fReferenceSolution = fSolution;
	Assemble();
    TPZAbstractFrontMatrix<STATE> *solverMat = dynamic_cast<TPZAbstractFrontMatrix<STATE> *>(MatrixSolver<STATE>().Matrix().operator->());
	//    fSolver->Solve(fRhs, fRhs);
    if(solverMat->IsDecomposed() == ELU)
    {
        TPZFMatrix<STATE> &rhs = fRhs;
        solverMat->Subst_Forward(&rhs);
    } else if(solverMat->IsDecomposed() == ECholesky)
    {
        TPZFMatrix<STATE> &rhs = fRhs;
        solverMat->Subst_Forward(&rhs);
    } else if(solverMat->IsDecomposed() == ELDLt)
    {
        PZError << __PRETTY_FUNCTION__;
        PZError << ": dont know what to do...\n";
        DebugStop();
    }    
}
template<class TVar>
void TPZSubMeshFrontalAnalysis::CondensedSolution(TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	//	ek = fReducableStiff.K11Red();
	//	ef = fReducableStiff.F1Red();
	//ek.Print("ek condensed");
	if(fFront) {
        TPZFMatrix<TVar> &rhs = fRhs;
		fFront->ExtractFrontMatrix(ek);
		int64_t next = ek.Rows();
		int64_t neq = fRhs.Rows();
		ef.Redim(next,1);
		int eq;
		for(eq=0; eq<next; eq++) {
			ef(eq,0) = rhs(eq+neq-next,0);
		}
	}
}

template<class TVar>
void TPZSubMeshFrontalAnalysis::LoadSolutionInternal(
    TPZFMatrix<TVar> &mySol,const TPZFMatrix<TVar> &myRhs,
    const TPZFMatrix<TVar> &myRefSol, const TPZFMatrix<TVar> &sol)
{
    int64_t numinter = fMesh->NumInternalEquations();
	int64_t numeq = fMesh->TPZCompMesh::NEquations();
	TPZFMatrix<TVar> soltemp(numeq,1,0.);
    TPZAbstractFrontMatrix<STATE> *solverMat = dynamic_cast<TPZAbstractFrontMatrix<STATE> *>(MatrixSolver<STATE>().Matrix().operator->());
	int64_t i;
	for(i=0;i<numinter;i++) soltemp(i,0) = myRhs.GetVal(i,0);
	for(; i<numeq; i++) {
		soltemp(i,0) = sol.GetVal(i,0)-myRefSol.GetVal(i,0);
	}
    if(solverMat->IsDecomposed() == ELU)
    {
        solverMat->Subst_Backward(&soltemp);
    } else if(solverMat->IsDecomposed() == ECholesky)
    {
        solverMat->Subst_Backward(&soltemp);
    } else if(solverMat->IsDecomposed() == ELDLt)
    {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" logic error. Aborting...\n";
        DebugStop();
    }
	
	mySol = myRefSol + soltemp;
	TPZLinearAnalysis::LoadSolution();
}

void TPZSubMeshFrontalAnalysis::LoadSolution(const TPZFMatrix<STATE> &sol)
{
	LoadSolutionInternal<STATE>(fSolution,fRhs,fReferenceSolution,sol);
}


template
void TPZSubMeshFrontalAnalysis::CondensedSolution<STATE>(TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
//TODOCOMPLEX
//template
// void TPZSubMeshFrontalAnalysis::CondensedSolution<CSTATE>(TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
