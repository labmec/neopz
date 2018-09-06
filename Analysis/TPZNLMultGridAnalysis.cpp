/**
 * @file
 * @brief Contains implementations of the TPZNonLinMultGridAnalysis methods.
 */
#include "TPZNLMultGridAnalysis.h"
#include <math.h>                // for fabs
#include <stdio.h>               // for NULL
#include <string.h>              // for strcmp
#include <time.h>                // for clock, CLOCKS_PER_SEC
#include <iostream>              // for operator<<, string, basic_ostream, cout
#include <map>                   // for map, map<>::iterator, __map_iterator
#include <utility>               // for pair
#include "TPZCompElDisc.h"       // for TPZCompElDisc
#include "pzadmchunk.h"          // for TPZAdmChunkVector
#include "TPZChunkVector.h"             // for TPZChunkVector
#include "pzcmesh.h"             // for TPZCompMesh
#include "pzcompel.h"            // for TPZCompEl
#include "pzconslaw.h"           // for TPZConservationLaw
#include "pzdxmesh.h"            // for TPZDXGraphMesh
#include "pzvtkmesh.h"            // for TPZDXGraphMesh
#include "pzeltype.h"            // for MElementType::EDiscontinuous, MEleme...
#include "pzerror.h"             // for PZError
#include "pzflowcmesh.h"         // for TPZFlowCompMesh
#include "pzfunction.h"          // for TPZFunction
#include "pzgeoel.h"             // for TPZGeoEl
#include "pzgmesh.h"             // for TPZGeoMesh
#include "TPZMaterial.h"          // for TPZMaterial
#include "pzskylstrmatrix.h"     // for TPZSkylineStructMatrix
#include "pzsolve.h"             // for TPZMatrixSolver, TPZSolver, TPZMatri...
#include "pzstepsolver.h"        // for TPZStepSolver
#include "pzstrmatrix.h"         // for TPZStructMatrixOR
#include "pztransfer.h"          // for TPZTransfer
#include "pzvec.h"               // for TPZVec
#include "tpzagglomeratemesh.h"  // for TPZAgglomerateMesh
#include "tpzautopointer.h"      // for TPZAutoPointer
#include "TPZAgglomerateEl.h"    // for TPZAgglomerateElement
using namespace std;


TPZCompMesh *TPZNonLinMultGridAnalysis::IMesh(int64_t index){
	
	if( index < 0 || index >= fMeshes.NElements() )
		PZError << "TPZNonLinMultGridAnalysis::IMesh mesh index out of range\n";
	return fMeshes[index];
}

TPZNonLinMultGridAnalysis::TPZNonLinMultGridAnalysis(TPZCompMesh *cmesh) : 
TPZAnalysis(cmesh), fBegin(0), fInit(0) {
	cmesh->SetName("* * * MALHA INICIAL * * *");
	fMeshes.Push(cmesh);
	TPZStepSolver<STATE> solver;
	solver.SetDirect(ELDLt);
	TPZMatrixSolver<STATE> *clone = dynamic_cast<TPZMatrixSolver<STATE> *>(solver.Clone());
	SetSolver(*clone);
	fSolvers.Push(clone);
	fSolutions.Push(new TPZFMatrix<STATE>(fSolution));
	fPrecondition.Push(0);
}

TPZNonLinMultGridAnalysis::~TPZNonLinMultGridAnalysis() {
	while (fMeshes.NElements()) delete fMeshes.Pop();
	while(fSolutions.NElements()) delete fSolutions.Pop();
	while(fSolvers.NElements()) delete fSolvers.Pop();
	while(fPrecondition.NElements()) delete fPrecondition.Pop();
}


void TPZNonLinMultGridAnalysis::AppendMesh(TPZCompMesh * mesh){
	
	if(fMeshes.NElements() != fSolvers.NElements() || fMeshes.NElements() != fSolutions.NElements() ||
	   fPrecondition.NElements() != fMeshes.NElements()) {
		cout << "TPZNonLinMultGridAnalysis::AppendMesh can only be called after solving the coarse mesh\n";
		return;
	}
	fMeshes.Push(mesh);
}

TPZCompMesh *TPZNonLinMultGridAnalysis::PopMesh() {
	
	if(fMeshes.NElements() == 1) {
		cout << "TPZNonLinMultGridAnalysis cannot delete the root mesh, sorry\n";
		return 0;
	}
	if(fSolutions.NElements() == fMeshes.NElements()) delete fSolutions.Pop();
	delete fSolvers.Pop();
	SetSolver(*fSolvers[fSolvers.NElements()-1]);
	fCompMesh = fMeshes[fMeshes.NElements()-2];
	fSolution = *fSolutions[fSolutions.NElements()-1];
	return fMeshes.Pop();
}

TPZCompMesh *TPZNonLinMultGridAnalysis::AgglomerateMesh(TPZCompMesh *finemesh,
														int levelnumbertogroup){
	
	TPZVec<int64_t> accumlist;
	int dim = finemesh->Dimension();
	int64_t numaggl;
	TPZAgglomerateElement::ListOfGroupings(finemesh,accumlist,levelnumbertogroup,numaggl,dim);
	TPZCompMesh *aggmesh;
	aggmesh = TPZAgglomerateElement::CreateAgglomerateMesh(finemesh,accumlist,numaggl);
	return aggmesh;
}

TPZCompMesh  *TPZNonLinMultGridAnalysis::UniformlyRefineMesh(TPZCompMesh *coarcmesh,int levelnumbertorefine,int setdegree) {
	
	if(levelnumbertorefine < 1) return coarcmesh;
	TPZGeoMesh *gmesh = coarcmesh->Reference();
	if(!gmesh) {
		cout << "TPZMGAnalysis::UniformlyRefineMesh mesh with null reference, cancelled method\n";
		return 0;
	}
	cout << "\nTPZNonLinMultGridAnalysis::UniformlyRefineMesh uniforme division of coarcmesh,"
	<< " levels to be fine = " << levelnumbertorefine << endl;
	
	gmesh->ResetReference();
	TPZCompMesh *finemesh;
	
	finemesh = new TPZFlowCompMesh(gmesh);
	
	std::map<int, TPZMaterial * >::iterator m;
	for(m=coarcmesh->MaterialVec().begin(); m != coarcmesh->MaterialVec().end(); m++) {
		TPZMaterial * mat = m->second;
		if(!mat) continue;
		mat->Clone(finemesh->MaterialVec());
	}
	TPZAdmChunkVector<TPZCompEl *> &elementvec = coarcmesh->ElementVec();
	int64_t el,nelem = elementvec.NElements();
	for(el=0; el<nelem; el++) {
		TPZCompEl *cel = elementvec[el];
		if(!cel) continue;
		if(cel->Type() == EAgglomerate){
			PZError << "TPZNonLinMultGridAnalysis::UniformlyRefineMesh mesh error,"
			<< " not refined\n";
			return new TPZCompMesh(NULL);
		}
		if(cel->Type() != EDiscontinuous) continue;
		TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
		int degree = disc->Degree();
		
		TPZGeoEl *gel = disc->Reference();
		if(!gel) {
			cout << "TPZMGAnalysis::UniformlyRefineMesh encountered an element without"
			<< " geometric reference\n";
			continue;
		}
		TPZStack<TPZGeoEl *> sub,sub1;
		//GetRefinedGeoEls(geo,sub);
		int lev = 0,k,nsons,i;
		gel->Divide(sub);
		lev++;
		while(lev <  levelnumbertorefine){
			int nsubs = sub.NElements();
			TPZVec<TPZGeoEl *> copy(sub);
			sub.Resize(0);
			for(i=0;i<nsubs;i++){
				copy[i]->Divide(sub1);
				nsons = sub1.NElements();
				for(k=0;k<nsons;k++) sub.Push(sub1[k]);
			}
			lev++;
		}
		int nsub = sub.NElements(),isub;
		int64_t index;
		//o construtor adequado ja deveria ter sido definido
		for(isub=0; isub<nsub; isub++) {
			disc = dynamic_cast<TPZCompElDisc *>(finemesh->CreateCompEl(sub[isub],index));
			if(setdegree > 0 && setdegree != degree) disc->SetDegree(degree);
			//caso setdegree < 0 preserva-se o grau da malha inicial
		}
	}
	finemesh->AdjustBoundaryElements();
	return finemesh;
}

void TPZNonLinMultGridAnalysis::ResetReference(TPZCompMesh *aggcmesh){
	//APLICAR ESTA FUN�O ANTES DE GERAR A MALHA COM O DX
	
	//caso o aglomerado tem refer�cia anulam-se as referencias
	//dos sub-elementos 'geom�ricos' aglomerados por ele
	//caso contr�io deixa-se um nico elemento geom�rico
	//apontando para o aglomerado
	//isso forma uma parti�o da malha atual por elementos computacionais
	
	int64_t nel = aggcmesh->NElements(),i;
	TPZCompMesh *finemesh;
	//n� todo index �sub-elemento
	for(i=0;i<nel;i++){
		TPZCompEl *cel = aggcmesh->ElementVec()[i];
		if(!cel) continue;
		if(cel->Type() == EInterface) continue;
		if(cel->Type() == EDiscontinuous) continue;
		TPZMaterial * mat = cel->Material();
		if(!mat) PZError << "TPZNonLinMultGridAnalysis::ResetReference null material\n";
		if(mat->Id() < 0) continue;
		TPZGeoEl *gel = cel->Reference();
		TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
		if(!agg) 
			PZError << "TPZNonLinMultGridAnalysis::ResetReference not agglomerate element\n";
		finemesh = agg->MotherMesh();
		if(!finemesh) 
			PZError << "TPZNonLinMultGridAnalysis::ResetReference null fine mesh\n";
		TPZStack<int64_t> vec;
		agg->IndexesDiscSubEls(vec);
		int64_t i,size = vec.NElements();
		if(!size) PZError << "main::ResetReference error1\n";
		TPZCompEl *sub0 = finemesh->ElementVec()[vec[0]],*sub;
		for(i=1;i<size;i++){
			sub = finemesh->ElementVec()[vec[i]];
			TPZGeoEl *ref = sub->Reference();
			if(!ref) PZError << "main::ResetReference error2\n";
			ref->SetReference(NULL);
			//o aglomerado n� tem geom�rico direto associado
			//agora existe um nico geom�rico apontando
			//para ele
		}
		if(gel){
			TPZGeoEl * ref0 = sub0->Reference();
			if(!ref0) PZError << "main::ResetReference error2\n";
			ref0->SetReference(NULL);
			//o aglomerado tem geom�rico direto associado
			//e esse aponta para ele
		}
	}
}

void TPZNonLinMultGridAnalysis::SetReference(TPZCompMesh *aggcmesh){
	
	int64_t nel = aggcmesh->NElements(),i;
	//TPZCompMesh *finemesh;
	//n� todo index �sub-elemento
	for(i=0;i<nel;i++){
		TPZCompEl *cel = aggcmesh->ElementVec()[i];
		if(!cel) continue;
		if(cel->Type() == EInterface) continue;
		if(cel->Type() == EDiscontinuous) continue;
		TPZMaterial * mat = cel->Material();
		if(!mat) PZError << "TPZNonLinMultGridAnalysis::SetReference null material\n";
		if(mat->Id() < 0) continue;
		TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
		if(!agg) 
			PZError << "TPZNonLinMultGridAnalysis::SetReference not agglomerate element\n";
		TPZStack<int64_t> elvec;
		agg->IndexesDiscSubEls(elvec);
		//os computacionais da malha fina apontam para os respectivos geometricos
		//os geometricos deveram apontar para o agglomerado que o agrupa;
		//si existe um geometrico tal que as referencias dos agrupados no aglomerado
		//formam uma particao unitaria desse entao esse geometrico ja
		//aponta para esse aglomerado
		int64_t indsize = elvec.NElements(),k;
		for(k=0;k<indsize;k++){
			TPZCompEl *cel = agg->MotherMesh()->ElementVec()[elvec[k]];
			if(!cel){
				PZError << "TPZNonLinMultGridAnalysis::SetReference null sub-element\n";
				continue;
			}
			TPZGeoEl *ref = cel->Reference();
			ref->SetReference(agg);
		}
	}
}

void TPZNonLinMultGridAnalysis::CoutTime(clock_t &start,const char *title){
    clock_t end = clock();
    cout << title <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

void TPZNonLinMultGridAnalysis::SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial * mat){
	
	TPZFlowCompMesh *fm  = dynamic_cast<TPZFlowCompMesh *>(CompMesh);
	TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
	REAL timestep = law->TimeStep();
	if(timestep <= 0.0){
		fFunction(mat,CompMesh);
		return;
	}
	//int nstate = mat->NStateVariables();
	
	REAL maxveloc;//,gama = 1.4;
	//maxveloc = fm->MaxVelocityOfMesh(nstate,gama);
	maxveloc = fm->MaxVelocityOfMesh();
	REAL deltax = CompMesh->LesserEdgeOfMesh();
	//REAL deltax = CompMesh->DeltaX();
	//REAL deltax = CompMesh->MaximumRadiusOfEl();
	int degree = CompMesh->GetDefaultOrder();
	law->SetTimeStep(maxveloc,deltax,degree);   //JorgeC
}

void TPZNonLinMultGridAnalysis::SmoothingSolution(REAL tol,int numiter,TPZMaterial * mat,TPZAnalysis &an,TPZFMatrix<STATE> &rhs){
	
	//pelo menos duas iteracoes para calcular o residuo
	if(numiter <= 1) numiter = 2;
	TPZCompMesh *anmesh = an.Mesh();
	int iter = 0;
	cout << "PZAnalysis::SmoothingSolutionTest beginning of the iterative process iterac� = " << iter << endl;
	an.Solution().Zero();
	SetDeltaTime(anmesh,mat);
	an.Run();
	cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteration = " << ++iter << endl;
	an.LoadSolution();
	mat->SetForcingFunction(0);
	REAL normsol = Norm(Solution());
	TPZFMatrix<STATE> rhsim1 = an.Rhs();
	
	while(iter < numiter && normsol < tol) {
		
		an.Solution().Zero();
		SetDeltaTime(anmesh,mat);
		//    rhsim1 = an.Rhs();//guarda o anterior res�uo
		an.Run();
		an.LoadSolution();
		cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteracao = " << ++iter << endl;
		normsol = Norm(Solution());
	}
	if(iter < numiter){
		cout << "\nTPZNonLinMultGridAnalysis::SmoothingSolution the iterative process stopped"
		<< " due the great norm of the solution, norm solution = " << normsol << endl;
	}
	an.LoadSolution();
	rhs = an.Rhs() - rhsim1;//ltimo residuo - res�uo anterior
}


void TPZNonLinMultGridAnalysis::SmoothingSolution(REAL tol,int numiter,TPZMaterial * mat,
												  TPZAnalysis &an,int marcha,const std::string &dxout) {
	if(numiter <= 0){
		cout << "PZAnalysis::SmoothingSolution NENHUMA RESOLU�O EFETUADA\n";
		an.Solution().Zero();
		return;
	}
	if(marcha){
		//saida com o DX
		SmoothingSolution2(tol,numiter,mat,an,marcha,dxout);
		return;
	}
	TPZCompMesh *anmesh = an.Mesh();
	cout << "PZAnalysis::SmoothingSolutionTest beginning of the iterative process\n";
	int iter = 0;
	an.Solution().Zero();
	SetDeltaTime(anmesh,mat);
	an.Run();
	cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteration = " << ++iter << endl;
	an.LoadSolution();
	mat->SetForcingFunction(0);
	REAL normsol = Norm(Solution());
	
	while(iter < numiter && normsol < tol) {
		
		an.Solution().Zero();
		SetDeltaTime(anmesh,mat);
		an.Run();
		an.LoadSolution();
		cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteracao = " << ++iter << endl;
		normsol = Norm(Solution());
	}
	if(iter < numiter){
		cout << "\nTPZNonLinMultGridAnalysis::SmoothingSolution the iterative process stopped"
		<< " due the great norm of the solution, norm solution = " << normsol << endl;
	}
	an.LoadSolution();
}

void TPZNonLinMultGridAnalysis::SmoothingSolution2(REAL tol,int numiter,TPZMaterial * mat,
												   TPZAnalysis &an,int marcha,const std::string &dxout) {
	
	TPZVec<std::string> scalar(1),vector(0);
	scalar[0] = "pressure";
	int dim = mat->Dimension();
	TPZCompMesh *anmesh = an.Mesh();
//	ResetReference(anmesh);//retira refer�cias para criar graph consistente
	TPZVTKGraphMesh graph(anmesh, dim, mat, scalar, vector);
//	TPZDXGraphMesh graph(anmesh,dim,mat,scalar,vector);
//	SetReference(anmesh);//recupera as refer�cias retiradas
	graph.SetFileName(dxout);
	int resolution = 0;
	graph.SetResolution(resolution);
	graph.DrawMesh(dim);
	int iter = 0,draw=0;
	fInit = clock();
	an.Solution().Zero();
	fBegin = clock();
	SetDeltaTime(anmesh,mat);
	cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteration = " << ++iter
	<< " general time 0\n";
	an.Run();
	CoutTime(fBegin,"TPZNonLinMultGridAnalysis:: Fim system solution first iteration");
	fBegin = clock();
	an.LoadSolution();
	REAL time = (REAL)iter;
	graph.DrawSolution(draw++,time);
	fBegin = clock();
	mat->SetForcingFunction(0);
	REAL normsol = Norm(Solution());
	
	while(iter < numiter && normsol < tol) {
		
		fBegin = clock();
		an.Solution().Zero();
		SetDeltaTime(anmesh,mat);
		an.Run();
		CoutTime(fBegin,"TPZNonLinMultGridAnalysis:: Fim system solution actual iteration");
		CoutTime(fInit,"TPZNonLinMultGridAnalysis:: accumulated time");
		fBegin = clock();
		an.LoadSolution();
		cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteracao = " << ++iter << endl;
		normsol = Norm(Solution());
		if( REAL(iter) / REAL(marcha) == draw || marcha == 1){
			time = (REAL)iter;
			graph.DrawSolution(draw++,time);
		}
	}
	if(iter < numiter){
		cout << "\nTPZNonLinMultGridAnalysis::SmoothingSolution the iterative process stopped"
		<< " due the great norm of the solution, norm solution = " << normsol << endl;
	}
	an.LoadSolution();
	CoutTime(fInit,"TPZNonLinMultGridAnalysis::SmoothingSolution general time of iterative process");
	time = (REAL)iter;
	graph.DrawSolution(draw++,time);
	ofstream out("SmoothingSolution2_ANALYSIS.out");
	an.Print("\n\n* * * SOLUCAO PELO ANALYSIS: ltima solu�o  * * *\n\n",out);
	out.flush();
	out.close();
}

void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix<STATE> &sol,TPZAnalysis &an,
											 const std::string &decompose,TPZFMatrix<STATE> &res){
	
	TPZAutoPointer<TPZMatrix<STATE> > stiff = an.Solver().Matrix();
	ofstream out("CalcResidual_STIFF.out");
	stiff->Print("\n\n\t\t\t* * * MATRIZ DE RIGIDEZ * * *\n\n",out);
	int dim = stiff->Dim(),i,j;
	res.Redim(dim,1);
	//c�culo de stiff * solution
	TPZFMatrix<STATE> tsup(dim,1),diag(dim,1),tinf(dim,1);
	
	if( !strcmp(decompose.c_str() , "LDLt") ) {
		//tri�gulo superior
		for(i=0;i<dim;i++){
			STATE sum = 0.;
			for(j=i+1;j<dim;j++){
				sum += stiff->GetVal(i,j) * sol(j,0);
			}
			tsup(i,0) = sol(i,0) + sum;
		}
		//diagonal
		for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,0);
		//tri�gulo inferior
		for(i=0;i<dim;i++) {
			STATE sum = 0.;
			for(j=0;j<i;j++) {
				sum += stiff->GetVal(i,j) * diag(j,0);
			}
			res(i,0) = sum + diag(i,0);
		}
		
		if(0) {
			ofstream out("MATRIZES_DECOMPOSI�O.out");
			sol.Print("\n* * * sol * * *\n",out);
			tsup.Print("\n* * * tsup * * *\n",out);
			diag.Print("\n* * * diag * * *\n",out);
			tinf.Print("\n* * * tinf * * *\n",out);
			//rhs.Print("\n* * * rhs * * *\n",out);
			res.Print("\n* * * res * * *\n",out);
			stiff->Print("\n* * * stiff * * *\n",out);
			out.flush();
			out.close();
		}
		//cout << "\nNorma do res�uo: " << Norm(res) << endl;//ORIGINAL
	} else {
		cout << "TPZNonLinMultGridAnalysis::CalcResidual Calculation of the residue for this"
		<< " decomposition is not implemented, implements now!\n";
	}
}


void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix<STATE> &sol,TPZFMatrix<STATE> &anres,
											 TPZFMatrix<STATE> &res,TPZAnalysis &an,const std::string &decompose){
	
	TPZAutoPointer<TPZMatrix<STATE> > stiff = an.Solver().Matrix();
	ofstream out("CalcResidual_STIFF.out");
	stiff->Print("\n\n\t\t\t* * * MATRIZ DE RIGIDEZ * * *\n\n",out);
	int dim = stiff->Dim(),i,j;
	//c�culo de stiff * solution
	TPZFMatrix<STATE> tsup(dim,1),diag(dim,1),tinf(dim,1);
	
	if( !strcmp(decompose.c_str() , "LDLt") ){
		//tri�gulo superior
		for(i=0;i<dim;i++){
			STATE sum = 0.;
			for(j=i+1;j<dim;j++){
				sum += stiff->GetVal(i,j);
			}
			tsup(i,0) = sol(i,0) + sum;
		}
		//diagonal
		for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,0);
		//tri�gulo superior
		for(i=0;i<dim;i++){
			STATE sum = 0.;
			for(j=0;j<i;j++){
				sum += stiff->GetVal(i,j) * diag(i,0);
			}
			tinf(i,0) = sum + diag(i,0);
		}
		//diferenca (f - stiff * x)
		//TPZFMatrix<STATE> rhs = an.Rhs();
		for(i=0;i<dim;i++) res(i,0) = anres.GetVal(i,0) - tinf(i,0);
		sol.Print("\n* * * sol * * *\n",cout);
		anres.Print("\n* * * anres * * *\n",cout);
		tinf.Print("\n* * * tinf * * *\n",cout);
		res.Print("\n* * * residuo * * *\n",cout);
	} else {
		cout << "TPZNonLinMultGridAnalysis::CalcResidual Calculation of the residue for this"
		<< " decomposition is not implemented, implements now!\n";
	}
}




/////////////////////////////////////////////////////////////////////////////
////                                                                     ////
////               ALGORITMO SIMPLES COM UMA MALHA                       ////
////                                                                     ////
/////////////////////////////////////////////////////////////////////////////

void TPZNonLinMultGridAnalysis::OneGridAlgorithm(std::ostream &out,int nummat){
	//ALGORITMO SIMPLES A UMA MALHA
	int iter,marcha;
	cout << "TPZNonLinMultGridAnalysis::OneGridAlgorithm Name of the out dx OneGridAlgorithm.dx\n";
	int levelnumbertorefine = 1;
	cout << "TPZNonLinMultGridAnalysis:: nmero de n�eis a dividir: ";
	cin >> levelnumbertorefine;
	TPZCompMesh *coarcmesh = fMeshes[0];//malha grosseira inicial
	int setdegree = -1;//preserva o grau da malha inicial
	TPZCompMesh *finemesh = UniformlyRefineMesh(coarcmesh,levelnumbertorefine,setdegree);
	TPZMaterial * finemat = finemesh->FindMaterial(nummat);
	int meshdim = coarcmesh->Dimension();
	finemesh->SetDimModel(meshdim);
	finemesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL FINA * * *\n\n");
	finemesh->Reference()->ResetReference();
	finemesh->LoadReferences();
	TPZAnalysis finean(finemesh);
	TPZSkylineStructMatrix finestiff(finemesh);
	finean.SetStructuralMatrix(finestiff);
	TPZStepSolver<STATE> finesolver;
	finesolver.SetDirect(ELDLt);
	finean.SetSolver(finesolver);
	finean.Solution().Zero();
	SetDeltaTime(finemesh,finemat);//para calcular o passo e estimar o nmero de iterac�s
	TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(finemat);
	law->SetTimeStep(-1);//para obter o c�culo antes da primeira soluc�
	cout << "\nTPZNonLinMultGridAnalysis::OneGridAlgorithm Numero de iteracoes ? :\n";
	cin >> iter;
	cout << "\nTPZNonLinMultGridAnalysis::OneGridAlgorithm Marcha ? :\n";
	cin >> marcha;
	REAL sol_tol = 1.e1;//valor m�imo da ||solu�o||
	std::string solout("OneGridAlgorithm.vtk");
//	std::string solout("OneGridAlgorithm.dx");
	SmoothingSolution(sol_tol,iter,finemat,finean,marcha,solout);
}

/////////////////////////////////////////////////////////////////////////////
////                                                                     ////
////               ALGORITMO MULTIGRID A DUAS MALHAS                     ////
////                                                                     ////
/////////////////////////////////////////////////////////////////////////////

void TPZNonLinMultGridAnalysis::TwoGridAlgorithm(std::ostream &out,int nummat){
	
	ifstream INd("DADOS.in");
	TPZCompMesh *coarcmesh = fMeshes[0];//malha grosseira inicial
	TPZGeoMesh *geomesh = fMeshes[0]->Reference();//nica malha geom�rica
	int meshdim = coarcmesh->Dimension();
	//criando a malha fina
	int levelnumbertorefine = 1;
	cout << "TPZNonLinMultGridAnalysis:: nmero de n�eis a dividir: ";
	INd >> levelnumbertorefine;
	int setdegree = -1;//preserva o grau da malha inicial
	//newmesh = 0: coarcmesh se tornou a malha fina
	TPZCompMesh *finemesh = UniformlyRefineMesh(coarcmesh,levelnumbertorefine,setdegree);
	finemesh->SetDimModel(meshdim);
	finemesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL FINA * * *\n\n");
	//obtendo-se a malha menos fina por agrupamento
	int levelnumbertogroup = levelnumbertorefine - 1;//sera obtido por agrupamento o n�el 0
	cout << "TPZNonLinMultGridAnalysis:: nmero do n�el a ser agrupado: ";
	INd >> levelnumbertogroup;
	TPZCompMesh *aggmesh = AgglomerateMesh(finemesh,levelnumbertogroup);
	aggmesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL AGLOMERADA * * *\n\n");
	aggmesh->SetDimModel(meshdim);
	AppendMesh(aggmesh);
	aggmesh->Reference()->SetName("\n\t\t\t* * * MALHA GEOM�RICA REFINADA * * *\n\n");
	aggmesh->Reference()->Print(out);//malha geom�rica �uma s�
	//out << "\n\n\t\t\t* * * MALHA AGLOMERADA * * *\n\n";
	aggmesh->Print(out);
	//out << "\n\n\t\t\t* * * MALHA FINA * * *\n\n";
	finemesh->Print(out);
	out.flush();
	
	//analysis na malha aglomerada
	TPZAnalysis coarsean(fMeshes[1]);
	TPZSkylineStructMatrix coarsestiff(fMeshes[1]);
	coarsean.SetStructuralMatrix(coarsestiff);
	TPZStepSolver<STATE> coarsesolver;
	coarsesolver.SetDirect(ELDLt);
	TPZMatrixSolver<STATE> *clone = dynamic_cast<TPZMatrixSolver<STATE> *>(coarsesolver.Clone());
	coarsean.SetSolver(*clone);
	fSolvers.Push(clone);
	coarsean.Solution().Zero();
	fSolutions.Push(new TPZFMatrix<STATE>(coarsean.Solution()));
	fPrecondition.Push(0);
	//analysis na malha fina
	AppendMesh(finemesh);
	TPZAnalysis finean(fMeshes[2]);
	TPZSkylineStructMatrix finestiff(fMeshes[2]);
	finean.SetStructuralMatrix(finestiff);
	TPZStepSolver<STATE> finesolver;
	finesolver.SetDirect(ELDLt);
	clone = dynamic_cast<TPZMatrixSolver<STATE> *>(finesolver.Clone());
	finean.SetSolver(*clone);
	fSolvers.Push(clone);
	finean.Solution().Zero();
	fSolutions.Push(new TPZFMatrix<STATE>(finean.Solution()));
	fPrecondition.Push(0);
	//preparac� para aplicar m�odo multigrid a duas malhas
	int preiter,positer,premarcha,posmarcha;
	TPZMaterial * finemat = finemesh->FindMaterial(nummat);
	int finedim = finemat->Dimension();
	TPZMaterial * coarsemat = fMeshes[1]->FindMaterial(nummat);
	int coarsedim = coarsemat->Dimension();
	cout << "\nNumero de iteracoes pre-suavisamento? :\n";
	INd >> preiter;
	//preiter = 10;
	cout << "main:: Parametro marcha : \n";
	INd >> premarcha;
	//premarcha = 0;
	cout << "\nNumero de iteracoes p�-suavisamento? :\n";
	INd >> positer;
	//cout << "main:: Parametro marcha no p�-suavisamento :\n";
	//INd >> posmarcha;
	posmarcha = premarcha;
	fMeshes[2]->Reference()->ResetReference();
	fMeshes[2]->LoadReferences();
	//ofstream *dxout = new ofstream("PreSmoothingEuler.dx");//par�etro de SmootSolution
	//ofstream *dxout2 = new ofstream("PostSmoothingEuler.dx");
	// TRANSFER�CIA DE SOLU�ES
	TPZTransfer<STATE> transfer;
	fMeshes[2]->BuildTransferMatrixDesc(*fMeshes[1],transfer);
	TPZFMatrix<STATE> projectsol;
	REAL normsolfine = 0.0,normsolcoar = 1.e10,erro;
	REAL errsol = fabs(normsolcoar - normsolfine);
	REAL gridtol = 0.01;
	int64_t coarneq = fMeshes[1]->NEquations();
	int64_t fineneq = fMeshes[2]->NEquations();
	TPZFMatrix<STATE> finesol(fineneq,1,0.),fineres(fineneq,1),finesol0,projfinesol;
	TPZFMatrix<STATE> coarsesol(coarneq,1,0.),projfineres(coarneq,1),rhs(coarneq,1),frhsk;
	TPZFMatrix<STATE> finesolkeep,coarsesolkeep;
	int mgmaxiter = 100,mgiter = 0;
	TPZVec<std::string> scalar(1),vector(0);
	scalar[0] = "pressure";
	geomesh->ResetReference();
	fMeshes[2]->LoadReferences();
	TPZDXGraphMesh finegraph(fMeshes[2],finedim,finemat,scalar,vector);
	finegraph.SetFileName("FineSmoothing.dx");
	finegraph.SetResolution(0);
	finegraph.DrawMesh(finedim);
	geomesh->ResetReference();
	fMeshes[1]->LoadReferences();
	ResetReference(fMeshes[1]);//retira refer�cias para criar coarsegraph consistente  
	TPZDXGraphMesh coarsegraph(fMeshes[1],coarsedim,coarsemat,scalar,vector);
	SetReference(fMeshes[1]);//recupera as refer�cias
	coarsegraph.SetFileName("CoarseSmoothing.dx");
	coarsegraph.SetResolution(0);
	coarsegraph.DrawMesh(coarsedim);
	REAL sol_tol = 1.e15;//valor m�imo da ||solu�o||
	REAL time,skm1 = 1.0;
	cout << "TwoGridAlgorithm entre sk-1 : ";
	INd >> skm1;
	REAL skm1inv = 1.0/skm1;
	int draw = 0,residuo;
	cout << "TwoGridAlgorithm nmero m�imo de itera�es : ";
	INd >> mgmaxiter;
	cout << "TwoGridAlgorithm res�uo [1:Lk-1(u0k-1)] ou res�uo [2:fk-1] ? : ";
	INd >> residuo;
	fInit = clock();
	cout << "PZAnalysis::SmoothingSolutionTest beginning of the iterative process, time = 0\n";
	
	while( errsol > gridtol && mgiter < mgmaxiter){
		
		cout << "\nTwoGridAlgorithm: iterac� nmero = " << mgiter << "\n\n";
		geomesh->ResetReference();//por no final
		fMeshes[2]->LoadReferences();//por no final
		finesolkeep = finesol;//= 0 ou guarda a solu�o a seguir da itera�o anterior
		SmoothingSolution(sol_tol,preiter,finemat,finean,rhs);// PASSO 1
		finesol = finean.Solution();
		{
			time = (REAL)mgiter;//pode ser a soma dos passos de tempo
			if( REAL(mgiter) / REAL(premarcha) == draw || premarcha == 1 || (mgiter+1) == mgmaxiter){
				finegraph.DrawSolution(draw,time);
			}
			erro = Norm(finesol - finesolkeep);
			cout << "TwoGridAlgorithm: ||finesol(i) - finesol(i-1)|| = " << erro << endl;
		}
		transfer.TransferResidual(rhs,projfineres);// PASSO 2
		//ser�somado �residuo da malha grossa proveniente 
		//da solu�o projetada projfinesol
		projfineres *= skm1;
		fMeshes[1]->ProjectSolution(projfinesol);// PASSO 3
		fMeshes[1]->LoadSolution(projfinesol);
		geomesh->ResetReference();//istas duas linhas s� feitas
		fMeshes[1]->LoadReferences();//em ProjectSolution acima
		if(positer){
			//a solu�o projetada avan� na malha grossa
			SmoothingSolution(sol_tol,positer,coarsemat,coarsean);// PASSO 5
			projfinesol = coarsean.Solution();
		}
		//duas formas de calcular o pr�imo passo
		if(residuo == 1){// forma 1
			coarsean.Solution().Zero();
			fSolvers[1]->SetMatrix(0);
			fSolvers[1]->SetMatrix(coarsean.StructMatrix()->CreateAssemble(coarsean.Rhs(),NULL));
			CalcResidual(projfinesol,coarsean,"LDLt",rhs);//rhs = Stiffcoar * projfinesol
			//no assemble ser�adicionado +fRhs por isso: rhs = rhs - fRhs
			rhs = rhs + projfineres;
			coarsean.Rhs() = rhs;//res�uo
			coarsean.Solve();
		} else if(residuo == 2){// forma 2
			coarsean.Solution().Zero();
			fSolvers[1]->SetMatrix(0);
			//argumento 1 espera? adicionar o res�uo anterior ao atual calculado
			fSolvers[1]->SetMatrix(coarsean.StructMatrix()->CreateAssemble(coarsean.Rhs(),NULL));
			TPZFMatrix<STATE> coarres = coarsean.Rhs() + projfineres;
			coarsean.Rhs() = coarres;//res�uo
			coarsean.Solve();
		}
		coarsesolkeep = coarsesol;// = 0 ou solu�o grossa da itera�o anterior
		coarsesol = coarsean.Solution();//solu�o para o residuo na malha grossa
		transfer.TransferSolution(coarsesol-projfinesol,finesol0);
		finesol0 *= skm1inv;
		finesol0 = finesol + finesol0;
		fMeshes[2]->LoadSolution(finesol0);// PASSO 6
		{
			if( REAL(mgiter) / REAL(posmarcha) == draw || posmarcha == 1 || (mgiter+1) == mgmaxiter){
				coarsegraph.DrawSolution(draw,time);
				draw++;// <- deve ser: posmarcha = premarcha
			}
			erro = Norm(coarsesol - coarsesolkeep);
			cout << "||coarsesol(i) - coarsesol(i-1)|| = " << erro << endl;
		}
		mgiter++;
	}
	CoutTime(fInit,"TPZNonLinMultGridAnalysis::SmoothingSolution general time of iterative process");
}

