/**
 * \file
 * @brief Contains implementations of the TPZNonLinMultGridAnalysis methods.
 */

#include "TPZNLMultGridAnalysis.h"
#include "TPZCompElDisc.h"
#include "pzflowcmesh.h"
#include "TPZAgglomerateEl.h"
#include "pzflowcmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pztransfer.h"
#include "pzadmchunk.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontStructMatrix.h"
#include "pzmgsolver.h"
#include "pzseqsolver.h"
#include "pzstepsolver.h"
#include "pzquad.h"
#include "pzmaterial.h"
//#include "TPZConservationLaw.h"
#include "TPZDiffusionConsLaw.h"
#include "pzdxmesh.h"
#include "pzsolve.h"
//#include "pztempmat.h"
#include "tpzagglomeratemesh.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;
//class TPZTransfer;

TPZCompMesh *TPZNonLinMultGridAnalysis::IMesh(int index){
	
	if( index < 0 || index >= fMeshes.NElements() )
		PZError << "TPZNonLinMultGridAnalysis::IMesh mesh index out of range\n";
	return fMeshes[index];
}

TPZNonLinMultGridAnalysis::TPZNonLinMultGridAnalysis(TPZCompMesh *cmesh) : 
TPZAnalysis(cmesh), fBegin(0), fInit(0) {
	cmesh->SetName("* * * MALHA INICIAL * * *");
	fMeshes.Push(cmesh);
	TPZStepSolver<REAL> solver;
	solver.SetDirect(ELDLt);
	TPZMatrixSolver<REAL> *clone = dynamic_cast<TPZMatrixSolver<REAL> *>(solver.Clone());
	SetSolver(*clone);
	fSolvers.Push(clone);
	fSolutions.Push(new TPZFMatrix<REAL>(fSolution));
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
	
	TPZVec<int> accumlist;
	int numaggl,dim = finemesh->Dimension();
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
	
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator m;
	for(m=coarcmesh->MaterialVec().begin(); m != coarcmesh->MaterialVec().end(); m++) {
		TPZAutoPointer<TPZMaterial> mat = m->second;
		if(!mat) continue;
		mat->Clone(finemesh->MaterialVec());
	}
	TPZAdmChunkVector<TPZCompEl *> &elementvec = coarcmesh->ElementVec();
	int el,nelem = elementvec.NElements();
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
		int nsub = sub.NElements(),isub,index;
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
	
	int nel = aggcmesh->NElements(),i;
	TPZCompMesh *finemesh;
	//n� todo index �sub-elemento
	for(i=0;i<nel;i++){
		TPZCompEl *cel = aggcmesh->ElementVec()[i];
		if(!cel) continue;
		if(cel->Type() == EInterface) continue;
		if(cel->Type() == EDiscontinuous) continue;
		TPZAutoPointer<TPZMaterial> mat = cel->Material();
		if(!mat) PZError << "TPZNonLinMultGridAnalysis::ResetReference null material\n";
		if(mat->Id() < 0) continue;
		TPZGeoEl *gel = cel->Reference();
		TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
		if(!agg) 
			PZError << "TPZNonLinMultGridAnalysis::ResetReference not agglomerate element\n";
		finemesh = agg->MotherMesh();
		if(!finemesh) 
			PZError << "TPZNonLinMultGridAnalysis::ResetReference null fine mesh\n";
		TPZStack<int> vec;
		agg->IndexesDiscSubEls(vec);
		int i,size = vec.NElements();
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
	
	int nel = aggcmesh->NElements(),i;
	//TPZCompMesh *finemesh;
	//n� todo index �sub-elemento
	for(i=0;i<nel;i++){
		TPZCompEl *cel = aggcmesh->ElementVec()[i];
		if(!cel) continue;
		if(cel->Type() == EInterface) continue;
		if(cel->Type() == EDiscontinuous) continue;
		TPZAutoPointer<TPZMaterial> mat = cel->Material();
		if(!mat) PZError << "TPZNonLinMultGridAnalysis::SetReference null material\n";
		if(mat->Id() < 0) continue;
		TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
		if(!agg) 
			PZError << "TPZNonLinMultGridAnalysis::SetReference not agglomerate element\n";
		TPZStack<int> elvec;
		agg->IndexesDiscSubEls(elvec);
		//os computacionais da malha fina apontam para os respectivos geometricos
		//os geometricos deveram apontar para o agglomerado que o agrupa;
		//si existe um geometrico tal que as referencias dos agrupados no aglomerado
		//formam uma particao unitaria desse entao esse geometrico ja
		//aponta para esse aglomerado
		int indsize = elvec.NElements(),k;
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

void TPZNonLinMultGridAnalysis::SetDeltaTime(TPZCompMesh *CompMesh,TPZAutoPointer<TPZMaterial> mat){
	
	TPZFlowCompMesh *fm  = dynamic_cast<TPZFlowCompMesh *>(CompMesh);
	TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat.operator->());
	REAL timestep = law->TimeStep();
	if(timestep <= 0.0){
		fFunction(mat.operator->(),CompMesh);
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

void TPZNonLinMultGridAnalysis::SmoothingSolution(REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,TPZAnalysis &an,TPZFMatrix<REAL> &rhs){
	
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
	TPZFMatrix<REAL> rhsim1 = an.Rhs();
	
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


void TPZNonLinMultGridAnalysis::SmoothingSolution(REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,
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
	//  static int file = -1;
	//  file++;
	//ofstream out("SmoothingSolution_ANALYSIS.out");
	//   ofstream *out[2];
	//   char *name;
	//   if(!file) name = "SmoothingSolution_ANALYSIS1.out";
	//   if( file) name = "SmoothingSolution_ANALYSIS2.out";
	//   out[file] = new ofstream(name);
	//   an.Print("\n\n* * * SOLUCAO PELO ANALYSIS: primeira solu�o  * * *\n\n",*out[file]);
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
	//   an.Print("\n\n* * * SOLUCAO PELO ANALYSIS: ltima solu�o  * * *\n\n",*out[file]);
	//   out[file]->flush();
	//   out[file]->close();
}

void TPZNonLinMultGridAnalysis::SmoothingSolution2(REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,
												   TPZAnalysis &an,int marcha,const std::string &dxout) {
	
	TPZVec<std::string> scalar(1),vector(0);
	scalar[0] = "pressure";
	int dim = mat->Dimension();
	TPZCompMesh *anmesh = an.Mesh();
	ResetReference(anmesh);//retira refer�cias para criar graph consistente
	TPZDXGraphMesh graph(anmesh,dim,mat,scalar,vector);
	SetReference(anmesh);//recupera as refer�cias retiradas
	//ofstream dxout = new ofstream("SmoothingSolution2.dx");
	//cout << "\nTPZNonLinMultGridAnalysis::IterativeProcess out file : .dx\n";
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

void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix<REAL> &sol,TPZAnalysis &an,
											 const std::string &decompose,TPZFMatrix<REAL> &res){
	
	//TPZStepSolver *solver = dynamic_cast<TPZStepSolver *>(&an.Solver());
	//TPZMatrix<REAL> *stiff = solver->Matrix();
	TPZAutoPointer<TPZMatrix<REAL> > stiff = an.Solver().Matrix();
	ofstream out("CalcResidual_STIFF.out");
	stiff->Print("\n\n\t\t\t* * * MATRIZ DE RIGIDEZ * * *\n\n",out);
	int dim = stiff->Dim(),i,j;
	res.Redim(dim,1);
	//c�culo de stiff * solution
	TPZFMatrix<REAL> tsup(dim,1),diag(dim,1),tinf(dim,1);
	
	if( !strcmp(decompose.c_str() , "LDLt") ){
		//tri�gulo superior
		for(i=0;i<dim;i++){
			REAL sum = 0.;
			for(j=i+1;j<dim;j++){
				sum += stiff->GetVal(i,j) * sol(j,0);
			}
			tsup(i,0) = sol(i,0) + sum;
		}
		//diagonal
		for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,0);
		//tri�gulo inferior
		for(i=0;i<dim;i++){
			REAL sum = 0.;
			for(j=0;j<i;j++){
				sum += stiff->GetVal(i,j) * diag(j,0);
			}
			//tinf(i,0) = sum + diag(i,0);//ORIGINAL
			res(i,0) = sum + diag(i,0);
		}
		//diferenca (f - stiff * x)
		//TPZFMatrix<REAL> rhs = an.Rhs();//ORIGINAL
		//for(i=0;i<dim;i++) res(i,0) = rhs.GetVal(i,0) - tinf(i,0);//ORIGINAL
		
		if(0){
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


void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix<REAL> &sol,TPZFMatrix<REAL> &anres,
											 TPZFMatrix<REAL> &res,TPZAnalysis &an,const std::string &decompose){
	
	//TPZStepSolver *solver = dynamic_cast<TPZStepSolver *>(&an.Solver());
	//TPZMatrix<REAL> *stiff = solver->Matrix();
	TPZAutoPointer<TPZMatrix<REAL> > stiff = an.Solver().Matrix();
	ofstream out("CalcResidual_STIFF.out");
	stiff->Print("\n\n\t\t\t* * * MATRIZ DE RIGIDEZ * * *\n\n",out);
	int dim = stiff->Dim(),i,j;
	//c�culo de stiff * solution
	TPZFMatrix<REAL> tsup(dim,1),diag(dim,1),tinf(dim,1);
	
	if( !strcmp(decompose.c_str() , "LDLt") ){
		//tri�gulo superior
		for(i=0;i<dim;i++){
			REAL sum = 0.;
			for(j=i+1;j<dim;j++){
				sum += stiff->GetVal(i,j);
			}
			tsup(i,0) = sol(i,0) + sum;
		}
		//diagonal
		for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,0);
		//tri�gulo superior
		for(i=0;i<dim;i++){
			REAL sum = 0.;
			for(j=0;j<i;j++){
				sum += stiff->GetVal(i,j) * diag(i,0);
			}
			tinf(i,0) = sum + diag(i,0);
		}
		//diferenca (f - stiff * x)
		//TPZFMatrix<REAL> rhs = an.Rhs();
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


//   if(1){//TESTE 3
//     cout << "\n\nTESTE 3: Project Solution\n\n";
//     //ESTE TESTE MOSTRA QUE A SOLU�O �TRANSFERIDA DA 
//     //MALHA FINA PARA A MALHA GROSSA DE FORMA CONSISTENTE
//     //PARA QUALQUER GRAU DE INTERPOLA�O 
//     neq = fMeshes[2]->NEquations();
//     REAL c = sqrt(2.0);
//     TPZFMatrix<REAL> finesol(neq,1,c),finesol2(neq,1,0.);;//para grau 1
//     if(1){
//     //com grau 1 esta solu�o com o quadril�ero
//     int l;//[0,4]x[0,4] dividido em 4 elementos
//     for(l=0; l<4 ;l++) finesol(l,0) = 2.0;//como malha fina
//     for(l=12;l<16;l++) finesol(l,0) = 4.0;
//     for(l=24;l<28;l++) finesol(l,0) = 6.0;
//     for(l=36;l<40;l++) finesol(l,0) = 4.0;
// //     for(l=0;l<12;l++) if(!(l%3)) finesol(l,0) = 2.0;//como malha fina
// //     for(l=12;l<24;l++) if(!(l%3)) finesol(l,0) = 4.0;//e sem dividir com
// //     for(l=24;l<36;l++) if(!(l%3)) finesol(l,0) = 6.0;//malha grossa da solu��
// //     for(l=36;l<48;l++) if(!(l%3)) finesol(l,0) = 4.0;//exata  X+Y
//     }
//     if(0){//para grau 2
//     //com grau 2 esta solu�o com o quadril�ero
//     int l;//[0,4]x[0,4] dividido em 4 elementos
//     c = sqrt(2.0);
//     REAL c2 = c*c;
//     for(l=0;l<24;l++) {//como malha fina
//       if(l <  4)    {finesol(l,0) = 2.0;  continue;}
//       if(l <  8)    {finesol(l,0) = 2.*c; continue;}
//       if(l < 12)    {finesol(l,0) = c2;   continue;}
//       if(l < 16)    {finesol(l,0) = 2.*c; continue;}
//       if(l < 20)    {finesol(l,0) = 0.0;  continue;}
//       if(l < 24)    {finesol(l,0) = c2;   continue;}
//     }
//     for(l=24;l<48;l++) {
//       if(l < 28)    {finesol(l,0) = 10.;  continue;}
//       if(l < 32)    {finesol(l,0) = 2.*c; continue;}
//       if(l < 36)    {finesol(l,0) = c2;   continue;}
//       if(l < 40)    {finesol(l,0) = 6.*c; continue;}
//       if(l < 44)    {finesol(l,0) = 0.0;  continue;}
//       if(l < 48)    {finesol(l,0) = c2;   continue;}
//     }
//     for(l=48;l<72;l++) {
//       if(l < 52)    {finesol(l,0) = 18.;  continue;}
//       if(l < 56)    {finesol(l,0) = 6.*c; continue;}
//       if(l < 60)    {finesol(l,0) = c2;   continue;}
//       if(l < 64)    {finesol(l,0) = 6.*c; continue;}
//       if(l < 68)    {finesol(l,0) = 0.0;  continue;}
//       if(l < 72)    {finesol(l,0) = c2;   continue;}
//     }
//     for(l=72;l<96;l++) {
//       if(l < 76)    {finesol(l,0) = 10.;  continue;}
//       if(l < 80)    {finesol(l,0) = 6.*c; continue;}
//       if(l < 84)    {finesol(l,0) = c2;   continue;}
//       if(l < 88)    {finesol(l,0) = 2.*c; continue;}
//       if(l < 92)    {finesol(l,0) = 0.0;  continue;}
//       if(l < 96)    {finesol(l,0) = c2;   continue;}
//     }
//     }
//     fMeshes[2]->LoadSolution(finesol);
//     neq = fMeshes[1]->NEquations();
//     TPZFMatrix<REAL> coarsesol(neq,1,0.);
//     fMeshes[1]->ProjectSolution(coarsesol);
//     fMeshes[1]->LoadSolution(coarsesol);
//     finesol.Print("\n\n\t\t\t* * * finesol antes de projetada * * *\n\n",OUT);
//     coarsesol.Print("\n\n\t\t\t* * * coarsesol * * *\n\n",OUT);
//     transfer.TransferSolution(coarsesol,finesol2);
//     finesol.Print("\n\n\t\t\t* * * finesol depois de transferida * * *\n\n",OUT);
//     TPZFMatrix<REAL> residuo = finesol - finesol2;
//     residuo.Print("\n\n\t\t\t* * * RESIDUO: PROJETADA - TRANSFERIDA * * *\n\n",OUT);
//     int nel = fMeshes[1]->ElementVec().NElements(),i;
//     TPZVec<REAL> csi(3,0.),x(3);
//     TPZAgglomerateElement *aggel;
//     for(i=0;i<nel;i++){
//       TPZCompEl *cel = fMeshes[1]->ElementVec()[i];
//       if(!cel) continue;
//       if(cel->Type() != EAgglomerate) continue;
//       aggel = dynamic_cast<TPZAgglomerateElement *>(cel);
//       //basta 1 para teste
//       int size = aggel->NIndexes(),k;
//       for(k=0;k<size;k++){
// 	fMeshes[2]->Reference()->ResetReference();
// 	fMeshes[2]->LoadReferences();
// 	TPZCompEl *fineel0 = aggel->FineElement(k);
// 	TPZGeoEl *ref = fineel0->Reference();
// 	TPZManVector<REAL> sol;
// 	int var = 1;//solu�o u densidade
// 	fineel0->Solution(csi,var,sol);
// 	ref->X(csi,x);
// 	OUT << "Elemento -> " << ref->Id() << "\n";
// 	OUT << "x = " << x[0] << " " << x[1] << " " << x[2] << " : ";
// 	OUT << "solu�o fina -> " << sol[0] << "\n";
// 	fMeshes[1]->Reference()->ResetReference();
// 	fMeshes[1]->LoadReferences();
// 	aggel->Solution(x,-var,sol);//ponto real
// 	OUT << "\t\t\tsolu�o grossa -> " << sol[0] << "\n";}
//     }
//   }

//   if(0){//TESTE 2
//     cout << "\n\nTESTE 2: Transfer Solution\n\n";
//   //ESTE TESTE MOSTRA QUE A SOLU�O �TRANSFERIDA DA 
//   //MALHA GROSSA PARA A MALHA FINA DE FORMA EXATA
//   //PARA QUALQUER GRAU DE INTERPOLA�O
//   TPZFMatrix<REAL> coarsesol(neq,1,1.);
//   REAL val0 = 1.,val1 = 0.,val2 = 0.1,val3 = 1.0;//com CFL = 0
//   //cout << "TPZNonLinMultGridAnalysis::TwoGridAlgorithm: Entre valores : ";
//   //cin >> val1 >> val2 >> val3;
//   REAL eps = 0.001;
//   for(int i=0;i<neq;i++){
//     int k = (i+1)%4;
//     if(k == 2){
//       coarsesol(i,0) = val1;// + eps * i;
//       continue;
//     }
//     if(k == 3){
//       coarsesol(i,0) = val2 + eps * i;
//       continue;
//     }
//     if(k == 0) coarsesol(i,0) = val3 + eps * i;
//     if(k == 1) coarsesol(i,0) = val0 + eps * i;
//   }
//   neq = fMeshes[2]->NEquations();
//   TPZFMatrix<REAL> finesol(neq,1,0.);
//   fMeshes[1]->LoadSolution(coarsesol);
//   transfer.TransferSolution(coarsesol,finesol);
//   coarsesol.Print("\n\n\t\t\t* * * coarsesol * * *\n\n",OUT);
//   finesol.Print("\n\n\t\t\t* * * finesol * * *\n\n",OUT);
//   fMeshes[2]->LoadSolution(finesol);
//   int nel = fMeshes[1]->ElementVec().NElements(),i;
//   TPZVec<REAL> csi(3,0.),x(3);
//   TPZAgglomerateElement *aggel;
//   for(i=0;i<nel;i++){
//     TPZCompEl *cel = fMeshes[1]->ElementVec()[i];
//     if(!cel) continue;
//     if(cel->Type() != EAgglomerate) continue;
//     aggel = dynamic_cast<TPZAgglomerateElement *>(cel);
//     TPZStack<int> elvec;
//     aggel->IndexesDiscSubEls(elvec);
//     //basta 1 para teste
//     TPZCompEl *fineel0 = fMeshes[2]->ElementVec()[elvec[0]];
//     TPZGeoEl *ref = fineel0->Reference();
//     TPZManVector<REAL> sol;
//     int var = 4;//press�
//     fineel0->Solution(csi,var,sol);
//     ref->X(csi,x);
//     OUT << "Elemento -> " << ref->Id() << "\n";
//     OUT << "x = " << x[0] << " " << x[1] << " " << x[2] << " : ";
//     OUT << "solu�o fina -> " << sol[0] << "\n";
//     aggel->Solution(x,-var,sol);//ponto real
//     OUT << "\t\t\tsolu�o grossa -> " << sol[0] << "\n";
//   }
//   }

//   if(0){//TESTE 1: para verificar a qualidade da 
//         //         transfer�cia do res�uo
//   TPZFMatrix<REAL> coarsesol(neq,1,1.),transcoarsesol(neq,1,0.);
//   neq = fMeshes[2]->NEquations();
//   TPZFMatrix<REAL> finesol(neq,1,0.);
//   transfer.TransferSolution(coarsesol,finesol);
//   transfer.TransferResidual(finesol,transcoarsesol);
//   coarsesol.Print("\n\n\t\t\t* * * coarsesol * * *\n\n",OUT);
//   finesol.Print("\n\n\t\t\t* * * finesol * * *\n\n",OUT);
//   transcoarsesol.Print("\n\n\t\t\t* * * transcoarsesol * * *\n\n",OUT);
//   OUT.flush();
//   OUT.close();
//   }

//   if(0){
//   transfer.Print("\n\n\t\t\t* * * MATRIZ DE TRANSFER�CIA * * *\n\n",OUT);
//   int neq = fMeshes[1]->NEquations();
//   TPZFMatrix<REAL> coarseRhs(neq,1,0.);
//   transfer.TransferResidual(finean.Rhs(),coarseRhs);//transferindo r�iduo
//   coarseRhs.Print("\n\n* * * RES�UO TRANSFERIDO: FINO -> GROSSO * * *\n\n",OUT);
//   ofstream *dxout2 = new ofstream("PostSmoothingEuler.dx");
//   cout << "\nNumero de iteracoes p�-suavisamento? :\n";
//   cin >> numiter;
//   fMeshes[1]->Reference()->ResetReference();
//   fMeshes[1]->LoadReferences();
//   mat = fMeshes[1]->FindMaterial(nummat);
//   SmoothingSolution(tol,numiter,mat,coarsean,*dxout2,marcha);
//   coarsean.Rhs().Print("\n\n* * * RES�UO DA MALHA GROSSA * * *\n\n",OUT);
//   coarsean.Print("\n\n\t\t\t* * * ANALYSIS MALHA GROSSA ap� SmoothingSolution * * *\n\n",out);
//   neq = fMeshes[2]->NEquations();
//   TPZFMatrix<REAL> finesol(neq,1,0.);
//   transfer.TransferSolution(coarsean.Solution(),finesol);//transferindo solu�o
//   fMeshes[2]->LoadSolution(finesol);
//   TPZFMatrix<REAL> res;
//   CalcResidual(finesol,finean,"LDLt",res);
//   res.Print("\n\n\t\t\t* * * RES�UO MALHA FINA * * *\n\n",out);
//   finesol.Print("\n\n* * * SOLU�O TRANSFERIDA: GROSSA -> FINA * * *\n\n",OUT);

//   OUT.flush();
//   OUT.close();
//   }
// }

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
	TPZAutoPointer<TPZMaterial> finemat = finemesh->FindMaterial(nummat);
	int meshdim = coarcmesh->Dimension();
	finemesh->SetDimModel(meshdim);
	finemesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL FINA * * *\n\n");
	finemesh->Reference()->ResetReference();
	finemesh->LoadReferences();
	TPZAnalysis finean(finemesh);
	TPZSkylineStructMatrix finestiff(finemesh);
	finean.SetStructuralMatrix(finestiff);
	TPZStepSolver<REAL> finesolver;
	finesolver.SetDirect(ELDLt);
	finean.SetSolver(finesolver);
	finean.Solution().Zero();
	SetDeltaTime(finemesh,finemat.operator->());//para calcular o passo e estimar o nmero de iterac�s
	TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(finemat.operator->());
	law->SetTimeStep(-1);//para obter o c�culo antes da primeira soluc�
	cout << "\nTPZNonLinMultGridAnalysis::OneGridAlgorithm Numero de iteracoes ? :\n";
	cin >> iter;
	cout << "\nTPZNonLinMultGridAnalysis::OneGridAlgorithm Marcha ? :\n";
	cin >> marcha;
	REAL sol_tol = 1.e15;//valor m�imo da ||solu�o||
	std::string dxout("OneGridAlgorithm.dx");
	SmoothingSolution(sol_tol,iter,finemat,finean,marcha,dxout);
}

/////////////////////////////////////////////////////////////////////////////
////                                                                     ////
////               ALGORITMO MULTIGRID A DUAS MALHAS                     ////
////                                                                     ////
/////////////////////////////////////////////////////////////////////////////

void TPZNonLinMultGridAnalysis::TwoGridAlgorithm(std::ostream &out,int nummat){
	
	ifstream IN("DADOS.in");
	TPZCompMesh *coarcmesh = fMeshes[0];//malha grosseira inicial
	TPZGeoMesh *geomesh = fMeshes[0]->Reference();//nica malha geom�rica
	int meshdim = coarcmesh->Dimension();
	//criando a malha fina
	int levelnumbertorefine = 1;
	cout << "TPZNonLinMultGridAnalysis:: nmero de n�eis a dividir: ";
	IN >> levelnumbertorefine;
	int setdegree = -1;//preserva o grau da malha inicial
	//newmesh = 0: coarcmesh se tornou a malha fina
	TPZCompMesh *finemesh = UniformlyRefineMesh(coarcmesh,levelnumbertorefine,setdegree);
	finemesh->SetDimModel(meshdim);
	finemesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL FINA * * *\n\n");
	//obtendo-se a malha menos fina por agrupamento
	int levelnumbertogroup = levelnumbertorefine - 1;//sera obtido por agrupamento o n�el 0
	cout << "TPZNonLinMultGridAnalysis:: nmero do n�el a ser agrupado: ";
	IN >> levelnumbertogroup;
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
	TPZStepSolver<REAL> coarsesolver;
	coarsesolver.SetDirect(ELDLt);
	TPZMatrixSolver<REAL> *clone = dynamic_cast<TPZMatrixSolver<REAL> *>(coarsesolver.Clone());
	coarsean.SetSolver(*clone);
	fSolvers.Push(clone);
	coarsean.Solution().Zero();
	fSolutions.Push(new TPZFMatrix<REAL>(coarsean.Solution()));
	fPrecondition.Push(0);
	//analysis na malha fina
	AppendMesh(finemesh);
	TPZAnalysis finean(fMeshes[2]);
	TPZSkylineStructMatrix finestiff(fMeshes[2]);
	finean.SetStructuralMatrix(finestiff);
	TPZStepSolver<REAL> finesolver;
	finesolver.SetDirect(ELDLt);
	clone = dynamic_cast<TPZMatrixSolver<REAL> *>(finesolver.Clone());
	finean.SetSolver(*clone);
	fSolvers.Push(clone);
	finean.Solution().Zero();
	fSolutions.Push(new TPZFMatrix<REAL>(finean.Solution()));
	fPrecondition.Push(0);
	//preparac� para aplicar m�odo multigrid a duas malhas
	int preiter,positer,premarcha,posmarcha;
	TPZAutoPointer<TPZMaterial> finemat = finemesh->FindMaterial(nummat);
	int finedim = finemat->Dimension();
	TPZAutoPointer<TPZMaterial> coarsemat = fMeshes[1]->FindMaterial(nummat);
	int coarsedim = coarsemat->Dimension();
	cout << "\nNumero de iteracoes pre-suavisamento? :\n";
	IN >> preiter;
	//preiter = 10;
	cout << "main:: Parametro marcha : \n";
	IN >> premarcha;
	//premarcha = 0;
	cout << "\nNumero de iteracoes p�-suavisamento? :\n";
	IN >> positer;
	//cout << "main:: Parametro marcha no p�-suavisamento :\n";
	//IN >> posmarcha;
	posmarcha = premarcha;
	fMeshes[2]->Reference()->ResetReference();
	fMeshes[2]->LoadReferences();
	//ofstream *dxout = new ofstream("PreSmoothingEuler.dx");//par�etro de SmootSolution
	//ofstream *dxout2 = new ofstream("PostSmoothingEuler.dx");
	// TRANSFER�CIA DE SOLU�ES
	TPZTransfer transfer;
	fMeshes[2]->BuildTransferMatrixDesc(*fMeshes[1],transfer);
	TPZFMatrix<REAL> projectsol;
	REAL normsolfine = 0.0,normsolcoar = 1.e10,erro;
	REAL errsol = fabs(normsolcoar - normsolfine);
	REAL gridtol = 0.01;
	int coarneq = fMeshes[1]->NEquations();
	int fineneq = fMeshes[2]->NEquations();
	TPZFMatrix<REAL> finesol(fineneq,1,0.),fineres(fineneq,1),finesol0,projfinesol;
	TPZFMatrix<REAL> coarsesol(coarneq,1,0.),projfineres(coarneq,1),rhs(coarneq,1),frhsk;
	TPZFMatrix<REAL> finesolkeep,coarsesolkeep;
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
	IN >> skm1;
	REAL skm1inv = 1.0/skm1;
	int draw = 0,residuo;
	cout << "TwoGridAlgorithm nmero m�imo de itera�es : ";
	IN >> mgmaxiter;
	cout << "TwoGridAlgorithm res�uo [1:Lk-1(u0k-1)] ou res�uo [2:fk-1] ? : ";
	IN >> residuo;
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
			TPZFMatrix<REAL> coarres = coarsean.Rhs() + projfineres;
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

