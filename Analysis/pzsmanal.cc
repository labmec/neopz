// TPZSubMeshAnalysis.cpp: implementation of the TPZSubMeshAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#include "pzsmanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
#include "pztempmat.h"
#include "pzstrmatrix.h"
#include "pzsolve.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TPZSubMeshAnalysis::TPZSubMeshAnalysis(TPZSubCompMesh *mesh) : TPZAnalysis(mesh){
	fMesh = mesh;
}

TPZSubMeshAnalysis::~TPZSubMeshAnalysis()
{

}

void TPZSubMeshAnalysis::Assemble(){

	int numeq = fCompMesh->NEquations();
	int numinternal = fMesh->NumInternalEquations();
	fReferenceSolution.Redim(numeq,1);
	fRhs.Redim(numeq,1);
	fReducableStiff.Redim(numeq,numinternal);
	fSolver->SetMatrix(fStructMatrix->Create());
	fReducableStiff.SetK00(fSolver->Matrix());
	fReducableStiff.SetDecomposeType(ELDLt);
	fMesh->Assemble(fReducableStiff,fRhs);
	
}

void TPZSubMeshAnalysis::Run(ostream &out){

	//fReducableStiff.Print("Reducable stiff before assembled");
	fReferenceSolution = fSolution;
	Assemble();
	fReducableStiff.SetF(fRhs);
	//fReducableStiff.Print("Reducable stiff assembled");
	//fBlock->Print("Block structure",out);
	//fStiffness->Print("Stiffness matrix ",out);
	//out.flush();
	//fRhs->Print("Residual",out);
	//out.flush();
	
    //fStiffness->Print("Stiffness matrix ",out);
	//out.flush();
    //fRhs->Print("Residual",out);
    //out.flush();
}
void TPZSubMeshAnalysis::CondensedSolution(TPZFMatrix &ek, TPZFMatrix &ef){
	ek = fReducableStiff.K11Red();
	ef = fReducableStiff.F1Red();
	//ek.Print("ek condensed");
//	cout.flush();
}

void TPZSubMeshAnalysis::LoadSolution(TPZFMatrix &sol)
{

//	sol.Print("sol");
	int numinter = fMesh->NumInternalEquations();
	int numeq = fMesh->TPZCompMesh::NEquations();
	TPZFMatrix soltemp(numeq-numinter,1,0.);
//	fSolution->Print("Solucao a analise\n");
//	fReferenceSolution.Print("Solucao de referencia\n");
//	cout.flush();
	int i;
	for(i=0; i<numeq-numinter; i++) {
		soltemp(i,0) = sol(numinter+i,0)-fReferenceSolution(numinter+i,0);
	}
	TPZFMatrix uglobal(numeq,1,0.);
	fReducableStiff.UGlobal(soltemp,uglobal);
	fSolution = fReferenceSolution + uglobal;
//	fSolution->Print("fSolution");
	TPZAnalysis::LoadSolution();
//	cout.flush();
}
