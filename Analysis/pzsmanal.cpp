/**
 * \file
 * @brief Contains implementations of the TPZSubMeshAnalysis methods: implementation of the TPZSubMeshAnalysis class.
 */

#include "pzsmanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
//#include "pztempmat.h"
#include "pzstrmatrix.h"
#include "pzsolve.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TPZSubMeshAnalysis::TPZSubMeshAnalysis(TPZSubCompMesh *mesh) : TPZAnalysis(mesh){
	fMesh = mesh;
	fReducableStiff = new TPZMatRed<> ();
}

TPZSubMeshAnalysis::~TPZSubMeshAnalysis()
{
	
}

void TPZSubMeshAnalysis::Assemble(){
	
	std::cout << "Assembling the SubCompMesh index " << fMesh->Index() << std::endl;
	int numeq = fCompMesh->NEquations();
	int numinternal = fMesh->NumInternalEquations();
	fReferenceSolution.Redim(numeq,1);
	fRhs.Redim(numeq,1);
	fReducableStiff->Redim(numeq,numinternal);
	TPZMatRed<> *matred = dynamic_cast<TPZMatRed<> *> (fReducableStiff.operator->());
	fSolver->SetMatrix(fStructMatrix->Create());	
	//	fReducableStiff.SetK00(fSolver->Matrix());
	// this will initialize fK00 too
	matred->SetSolver(dynamic_cast<TPZMatrixSolver *>(fSolver->Clone()));
	//	TPZStructMatrix::Assemble(fReducableStiff,fRhs, *fMesh);
	time_t before = time (NULL);
	fStructMatrix->Assemble(fReducableStiff,fRhs,fGuiInterface);
	time_t after = time(NULL);
	double diff = difftime(after, before);
	std::cout << __PRETTY_FUNCTION__ << " tempo " << diff << std::endl;
}

void TPZSubMeshAnalysis::Run(std::ostream &out){
	
	//fReducableStiff.Print("Reducable stiff before assembled");
	fReferenceSolution = fSolution;
	time_t tempo = time(NULL);
	Assemble();
	time_t tempodepois = time(NULL);
	double elapsedtime = difftime(tempodepois, tempo);
	
	std::cout << "Tempo para assemblagem " << elapsedtime << std::endl;
	
	TPZMatRed<> *matred = dynamic_cast<TPZMatRed<> *> (fReducableStiff.operator->());
	matred->SetF(fRhs);
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
	time_t tempo = time(NULL);
	
	TPZMatRed<> *matred = dynamic_cast<TPZMatRed<> *> (fReducableStiff.operator->());
	ek = matred->K11Red();
	ef = matred->F1Red();
	//ek.Print("ek condensed");
	//	cout.flush();
	time_t tempodepois = time(NULL);
	double elapsedtime = difftime(tempodepois, tempo);
	
	std::cout << "Tempo para inversao " << elapsedtime << std::endl;
	
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
	TPZMatRed<> *matred = dynamic_cast<TPZMatRed<> *> (fReducableStiff.operator->());
	matred->UGlobal(soltemp,uglobal);
	fSolution = fReferenceSolution + uglobal;
	//	fSolution->Print("fSolution");
	TPZAnalysis::LoadSolution();
	//	cout.flush();
}
