/**
 * @file
 * @brief Contains implementations of the TPZSubMeshFrontalAnalysis methods.
 */

#include "pzsmfrontalanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
#include "TPZFrontStructMatrix.h"
#include "pzstrmatrix.h"
#include "pzsolve.h"

using namespace std;

// Construction/Destruction

TPZSubMeshFrontalAnalysis::TPZSubMeshFrontalAnalysis(TPZSubCompMesh *mesh) : TPZAnalysis(mesh)                                                , fReferenceSolution(false)//TODOCOMPLEX:set matrix type (complex/real)
{
	fMesh = mesh;
	fFront = 0;
}

TPZSubMeshFrontalAnalysis::~TPZSubMeshFrontalAnalysis()
{
	
}

void TPZSubMeshFrontalAnalysis::Run(std::ostream &out){
	
	//fReducableStiff.Print("Reducable stiff before assembled");
	fReferenceSolution = fSolution;
	Assemble();
	//    fSolver->Solve(fRhs, fRhs);
    if(fSolver->Matrix()->IsDecomposed() == ELU)
    {
        TPZFMatrix<STATE> &rhs = fRhs;
        fSolver->Matrix()->Subst_Forward(&rhs);
    } else if(fSolver->Matrix()->IsDecomposed() == ECholesky)
    {
        TPZFMatrix<STATE> &rhs = fRhs;
        fSolver->Matrix()->Subst_Forward(&rhs);
    } else if(fSolver->Matrix()->IsDecomposed() == ELDLt)
    {
        std::cout << "Dont know what to do...\n";
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
		int next = ek.Rows();
		int neq = fRhs.Rows();
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
    int numinter = fMesh->NumInternalEquations();
	int numeq = fMesh->TPZCompMesh::NEquations();
	TPZFMatrix<TVar> soltemp(numeq,1,0.);
	int i;
	for(i=0;i<numinter;i++) soltemp(i,0) = myRhs.GetVal(i,0);
	for(; i<numeq; i++) {
		soltemp(i,0) = sol.GetVal(i,0)-myRefSol.GetVal(i,0);
	}
    if(fSolver->Matrix()->IsDecomposed() == ELU)
    {
        fSolver->Matrix()->Subst_Backward(&soltemp);
    } else if(fSolver->Matrix()->IsDecomposed() == ECholesky)
    {
        fSolver->Matrix()->Subst_Backward(&soltemp);
    } else if(fSolver->Matrix()->IsDecomposed() == ELDLt)
    {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" logic error. Aborting...\n";
        DebugStop();
    }
	
	mySol = myRefSol + soltemp;
	TPZAnalysis::LoadSolution();
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