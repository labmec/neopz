// -*- c++ -*-

#include <stdlib.h>
#include "TPZNLMultGridAnalysis.h"
#include "TPZCompElDisc.h"
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
#include "TPZConservationLaw.h"
#include "pzonedref.h"
#include "pzdxmesh.h"
#include "pzsolve.h"
#include "pzflowcmesh.h"
using namespace std;
#include <string.h>
#include <stdio.h>
#include <iostream>
//class TPZTransfer;



TPZNonLinMultGridAnalysis::TPZNonLinMultGridAnalysis(TPZCompMesh *cmesh) : 
	TPZAnalysis(cmesh), fBegin(0), fInit(0) {
  fMeshes.Push(cmesh);
  TPZStepSolver solver;
  solver.SetDirect(ELDLt);
  SetSolver(solver);
  fSolutions.Push(&Solution());
  fSolvers.Push(&solver);
  fPrecondition.Push(&solver);
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
  int numaggl,dim = 2;
  TPZAgglomerateElement::ListOfGroupings(finemesh,accumlist,levelnumbertogroup,numaggl,dim);
  TPZCompMesh *aggmesh = new TPZFlowCompMesh(finemesh->Reference());
  TPZCompElDisc::CreateAgglomerateMesh(finemesh,*aggmesh,accumlist,numaggl);
  return aggmesh;
}


TPZCompMesh  *TPZNonLinMultGridAnalysis::UniformlyRefineMesh(TPZCompMesh *coarcmesh,int levelnumbertorefine,int setdegree) {

if(levelnumbertorefine < 1) return NULL;
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

  int nmat = coarcmesh->MaterialVec().NElements();
  int m;
  for(m=0; m<nmat; m++) {
    TPZMaterial *mat = coarcmesh->MaterialVec()[m];
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
      disc = dynamic_cast<TPZCompElDisc *>(sub[isub]->CreateCompEl(*finemesh,index));
      if(setdegree > 0 && setdegree != degree) disc->SetDegree(degree);
      //caso setdegree < 0 preserva-se o grau da malha inicial
    }
  }
  return finemesh;
}

void TPZNonLinMultGridAnalysis::ResetReference(TPZCompMesh *aggcmesh){
  //APLICAR ESTA FUNÇÃO ANTES DE GERAR A MALHA COM O DX

  //caso o aglomerado tem referência anulam-se as referencias
  //dos sub-elementos 'geométricos' aglomerados por ele
  //caso contrário deixa-se um único elemento geométrico
  //apontando para o aglomerado
  //isso forma uma partição da malha atual por elementos computacionais

  int nel = aggcmesh->NElements(),i;
  TPZCompMesh *finemesh;
  //não todo index é sub-elemento
  for(i=0;i<nel;i++){
    TPZCompEl *cel = aggcmesh->ElementVec()[i];
    if(!cel) continue;
    if(cel->Type() == EInterface) continue;
    if(cel->Type() == EDiscontinuous) continue;
    TPZMaterial *mat = cel->Material();
    if(!mat) PZError << "TPZIterativeAnalysis::ResetReference null material\n";
    if(mat->Id() < 0) continue;
    TPZGeoEl *gel = cel->Reference();
    TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
    if(!agg) PZError << "TPZIterativeAnalysis::ResetReference not agglomerate element\n";
    finemesh = agg->MotherMesh();
    if(!finemesh) PZError << "TPZIterativeAnalysis::ResetReference null fine mesh\n";
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
      //o aglomerado não tem geométrico direto associado
      //agora existe um único geométrico apontando
      //para ele
    }
    if(gel){
      TPZGeoEl * ref0 = sub0->Reference();
      if(!ref0) PZError << "main::ResetReference error2\n";
      ref0->SetReference(NULL);
      //o aglomerado tem geométrico direto associado
      //e esse aponta para ele
    }
  }
}

void TPZNonLinMultGridAnalysis::SetReference(TPZCompMesh *aggcmesh){

  int nel = aggcmesh->NElements(),i;
  //TPZCompMesh *finemesh;
  //não todo index é sub-elemento
  for(i=0;i<nel;i++){
    TPZCompEl *cel = aggcmesh->ElementVec()[i];
    if(!cel) continue;
    if(cel->Type() == EInterface) continue;
    if(cel->Type() == EDiscontinuous) continue;
    TPZMaterial *mat = cel->Material();
    if(!mat) PZError << "TPZIterativeAnalysis::SetReference null material\n";
    if(mat->Id() < 0) continue;
    TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
    if(!agg) PZError << "TPZIterativeAnalysis::SetReference not agglomerate element\n";
    TPZStack<int> elvec;
    agg->IndexesDiscSubEls(elvec);
    //os computacionais da malha fina apontam para os respectivos geométricos
    //os geométricos deveram apontar para o agglomerado que o agrupa;
    //si existe um geométrico tal que as referências dos agrupados no aglomerado
    //formam uma partição unitaria desse então esse geométrico já 
    //aponta para esse aglomerado
    int indsize = elvec.NElements(),k;
    for(k=0;k<indsize;k++){
      TPZCompEl *cel = agg->MotherMesh()->ElementVec()[elvec[k]];
      if(!cel){
	PZError << "TPZIterativeAnalysis::SetReference null sub-element\n";
	continue;
      }
      TPZGeoEl *ref = cel->Reference();
      ref->SetReference(agg);
    }
  }
}

void TPZNonLinMultGridAnalysis::CoutTime(clock_t &start,char *title){
    clock_t end = clock();
    cout << title <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

void TPZNonLinMultGridAnalysis::SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat){

  TPZFlowCompMesh *fm  = dynamic_cast<TPZFlowCompMesh *>(CompMesh);//= new TPZFlowCompMesh(CompMesh->Reference());
  REAL maxveloc = fm->MaxVelocityOfMesh();
  REAL deltax = CompMesh->LesserEdgeOfMesh();//REAL deltax = CompMesh->DeltaX();
  //REAL deltax = CompMesh->MaximumRadiusOfEl();
  TPZCompElDisc *disc;
  int degree = disc->gDegree;
  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
  law->SetDeltaTime(maxveloc,deltax,degree);
}

void TPZNonLinMultGridAnalysis::SmoothingSolution(REAL tol,int numiter,TPZMaterial *mat,TPZAnalysis &an) {

  fInit = clock();
  cout << "PZAnalysis::SmoothingSolutionTest beginning of the iterative process," 
       << " general time 0\n";

  //  int dim = mat->Dimension();
  int iter = 0;//,draw=0;
  an.Solution().Zero();
  fBegin = clock();
  an.Run();
  cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteration = " << ++iter << endl;
  CoutTime(fBegin,"TPZNonLinMultGridAnalysis:: Fim system solution first iteration");
  fBegin = clock();
  an.LoadSolution();
  SetDeltaTime(fCompMesh,mat);
  //  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
  fBegin = clock();
  mat->SetForcingFunction(0);
  REAL normsol = Norm(Solution());

  while(iter < numiter && normsol < tol) {
    
    fBegin = clock();
    an.Solution().Zero();
    an.Run();
    CoutTime(fBegin,"TPZNonLinMultGridAnalysis:: Fim system solution actual iteration");
    CoutTime(fInit,"TPZNonLinMultGridAnalysis:: accumulated time");
    fBegin = clock();
    an.LoadSolution();
    SetDeltaTime(fCompMesh,mat);
    cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteracao = " << ++iter << endl;
    normsol = Norm(Solution());
  }
  if(iter < numiter){
    cout << "\nTPZNonLinMultGridAnalysis::SmoothingSolution the iterative process stopped"
	 << " due the great norm of the solution, norm solution = " << normsol << endl;

  }
  an.LoadSolution();
  CoutTime(fInit,"TPZNonLinMultGridAnalysis::SmoothingSolution general time of iterative process");
}

void TPZNonLinMultGridAnalysis::TwoGridAlgorithm(ostream &out,int nummat){

  TPZCompMesh *coarcmesh = fMeshes[0];//malha grosseira inicial
  //criando a malha fina
  int levelnumbertorefine = 4;
  cout << "TPZNonLinMultGridAnalysis:: número de níveis a dividir: ";
  cin >> levelnumbertorefine;
  int setdegree = -1;//preserva o grau da malha inicial
  //newmesh = 0: coarcmesh se tornou a malha fina
  TPZCompMesh *finemesh = UniformlyRefineMesh(coarcmesh,levelnumbertorefine,setdegree);
  finemesh->SetDimModel(2);
  out << "\n\n\t\t* * * MALHA FINA * * *\n";
  finemesh->Reference()->Print(out);
  finemesh->Print(out);
  out.flush();
  //obtendo-se a malha menos fina por agrupamento
  int levelnumbertogroup = 2;//serão agrupados dois níveis de divisão
  cout << "TPZNonLinMultGridAnalysis:: número de níveis a agrupar: ";
  cin >> levelnumbertogroup;
  TPZCompMesh *aggmesh = AgglomerateMesh(finemesh,levelnumbertogroup);
  AppendMesh(aggmesh);
  out << "\n\n\t\t* * * MALHA AGLOMERADA * * *\n";
  aggmesh->Reference()->Print(out);
  aggmesh->Print(out);
  out.flush();
  //analysis na malha aglomerada
  TPZAnalysis coarsean(fMeshes[1]);
  TPZSkylineStructMatrix coarsestiff(fMeshes[1]);
  coarsean.SetStructuralMatrix(coarsestiff);
  coarsean.Solution().Zero();
  TPZStepSolver coarsesolver;
  coarsesolver.SetDirect(ELDLt);
  coarsean.SetSolver(coarsesolver);
  fSolutions.Push(&coarsean.Solution());
  fSolvers.Push(&coarsesolver);
  fPrecondition.Push(&coarsesolver);
  //analysis na malha fina
  AppendMesh(finemesh);
  TPZAnalysis finean(fMeshes[2]);
  TPZSkylineStructMatrix finestiff(fMeshes[2]);
  finean.SetStructuralMatrix(finestiff);
  finean.Solution().Zero();
  TPZStepSolver finesolver;
  finesolver.SetDirect(ELDLt);
  finean.SetSolver(finesolver);
  fSolutions.Push(&finean.Solution());
  fSolvers.Push(&finesolver);
  fPrecondition.Push(&finesolver);
  //suavisar a solução na malha fina
  REAL tol = 1.e15;
  int numiter = 100;
  int numeq = fMeshes[2]->NEquations();
  TPZFMatrix res(numeq,1);
  TPZMaterial *mat = fMeshes[2]->FindMaterial(nummat);
  SmoothingSolution(tol,numiter,mat,finean);
  //CalcResidual(finean.Solution(),res,finean,"LDLt");
  TPZTransfer transfer;
  finemesh->BuildTransferMatrixDesc(*aggmesh,transfer);
  SmoothingSolution(tol,numiter,mat,coarsean);
}

void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix &sol,TPZFMatrix &res,TPZAnalysis &an,char *decompose){

  //TPZStepSolver *solver = dynamic_cast<TPZStepSolver *>(&an.Solver());
  //TPZMatrix *stiff = solver->Matrix();
  TPZMatrix *stiff = an.Solver().Matrix();
  int dim = stiff->Dim(),i,j;
  //cálculo de stiff * solution
  TPZFMatrix tsup(dim,1),diag(dim,1),tinf(dim,1);

  if(decompose == "LDLt"){
    //triângulo superior
    for(i=0;i<dim;i++){
      REAL sum = 0.;
      for(j=i+1;j<dim;j++){
	 sum += stiff->GetVal(i,j);
      }
      tsup(i,0) = sol(i,0) + sum;
    }
    //diagonal
    for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,1);
    //triângulo superior
    for(i=0;i<dim;i++){
      REAL sum = 0.;
      for(j=0;j<i-1;j++){
	 sum += stiff->GetVal(i,j) * diag(i,0);
      }
      tinf(i,0) = sum + sol(i,0);
    }
    //diferenca (f - stiff * x)
    TPZMatrix rhs = an.Rhs();
    for(i=0;i<dim;i++) res(i,0) = Rhs()(i,0) - tinf(i,0);
  } else {
    cout << "TPZNonLinMultGridAnalysis::CalcResidual Calculation of the residue for this decomposition"
	 << " is not implemented, implements now!\n";
  }
}
