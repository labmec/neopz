/**
 * \file
 * @brief Contains implementations of the TPZSubMeshFrontalAnalysis methods.
 */

#include "pzsmfrontalanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
//#include "pztempmat.h"
#include "TPZFrontStructMatrix.h"
#include "pzstrmatrix.h"
#include "pzsolve.h"

using namespace std;

// Construction/Destruction

TPZSubMeshFrontalAnalysis::TPZSubMeshFrontalAnalysis(TPZSubCompMesh *mesh) : TPZAnalysis(mesh){
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
        fSolver->Matrix()->Subst_Forward(&fRhs);
    } else if(fSolver->Matrix()->IsDecomposed() == ECholesky)
    {
        fSolver->Matrix()->Subst_Forward(&fRhs);
    } else if(fSolver->Matrix()->IsDecomposed() == ELDLt)
    {
        std::cout << "Dont know what to do...\n";
        DebugStop();
    }
    
}
void TPZSubMeshFrontalAnalysis::CondensedSolution(TPZFMatrix &ek, TPZFMatrix &ef){
	//	ek = fReducableStiff.K11Red();
	//	ef = fReducableStiff.F1Red();
	//ek.Print("ek condensed");
	if(fFront) {
		fFront->ExtractFrontMatrix(ek);
		int next = ek.Rows();
		int neq = fRhs.Rows();
		ef.Redim(next,1);
		int eq;
		for(eq=0; eq<next; eq++) {
			ef(eq,0) = fRhs(eq+neq-next,0);
		}
	}
}

void TPZSubMeshFrontalAnalysis::LoadSolution(const TPZFMatrix &sol)
{
	
	//	sol.Print("sol");
	//	cout.flush();
	int numinter = fMesh->NumInternalEquations();
	int numeq = fMesh->TPZCompMesh::NEquations();
	TPZFMatrix soltemp(numeq,1,0.);
	//	fSolution.Print("Solucao a analise\n");
	//	fReferenceSolution.Print("Solucao de referencia\n");
	//	cout.flush();
	int i;
	for(i=0;i<numinter;i++) soltemp(i,0) = fRhs(i,0);
	for(; i<numeq; i++) {
		soltemp(i,0) = sol.GetVal(i,0)-fReferenceSolution(i,0);
	}
	//	soltemp.Print("Solucao temporaria");
	//	cout.flush();
    if(fSolver->Matrix()->IsDecomposed() == ELU)
    {
        fSolver->Matrix()->Subst_Backward(&soltemp);
    } else if(fSolver->Matrix()->IsDecomposed() == ECholesky)
    {
        fSolver->Matrix()->Subst_Backward(&soltemp);
    } else if(fSolver->Matrix()->IsDecomposed() == ELDLt)
    {
        std::cout << "Dont know what to do...\n";
        DebugStop();
    }
	
	fSolution = fReferenceSolution + soltemp;
	//	fSolution.Print("Final Solution");
	//	cout.flush();
	TPZAnalysis::LoadSolution();
}
