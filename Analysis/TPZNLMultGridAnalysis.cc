// -*- c++ -*-

#include "TPZNLMultGridAnalysis.h"
#include "TPZCompElDisc.h"
#include "TPZFlowCMesh.h"
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
#include "TPZDiffusionConsLaw.h"
#include "pzonedref.h"
#include "pzdxmesh.h"
#include "pzsolve.h"


#include "pztempmat.h"

using namespace std;
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
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
  TPZStepSolver solver;
  solver.SetDirect(ELDLt);
  TPZMatrixSolver *clone = dynamic_cast<TPZMatrixSolver *>(solver.Clone());
  SetSolver(*clone);
  fSolvers.Push(clone);
  fSolutions.Push(new TPZFMatrix(fSolution));
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
  TPZCompMesh *aggmesh = new TPZFlowCompMesh(finemesh->Reference());
  TPZCompElDisc::CreateAgglomerateMesh(finemesh,aggmesh,accumlist,numaggl);
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
  finemesh->AdjustBoundaryElements();
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
    if(!mat) PZError << "TPZNonLinMultGridAnalysis::SetReference null material\n";
    if(mat->Id() < 0) continue;
    TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
    if(!agg) 
      PZError << "TPZNonLinMultGridAnalysis::SetReference not agglomerate element\n";
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
	PZError << "TPZNonLinMultGridAnalysis::SetReference null sub-element\n";
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

void TPZNonLinMultGridAnalysis::SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat,int niter){

  TPZFlowCompMesh *fm  = dynamic_cast<TPZFlowCompMesh *>(CompMesh);
  if(!niter){
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
  TPZCompElDisc *disc;
  int degree = disc->gDegree;
  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
  law->SetDeltaTime(maxveloc,deltax,degree);
}

void TPZNonLinMultGridAnalysis::SmoothingSolution(REAL tol,int numiter,TPZMaterial *mat,
						  TPZAnalysis &an,ofstream &dxout,int marcha) {
  if(marcha){
    SmoothingSolution2(tol,numiter,mat,an,dxout,marcha);
    return;
  }
  TPZCompMesh *anmesh = an.Mesh();
  fInit = clock();
  cout << "PZAnalysis::SmoothingSolutionTest beginning of the iterative process," 
       << " general time 0\n";
  int iter = 0;
  an.Solution().Zero();
  fBegin = clock();
  SetDeltaTime(anmesh,mat,iter);
  an.Run();
  cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteration = " << ++iter << endl;
  CoutTime(fBegin,"TPZNonLinMultGridAnalysis:: Fim system solution first iteration");
  fBegin = clock();
  an.LoadSolution();
  static int file = -1;
  file++;
  //ofstream out("SmoothingSolution_ANALYSIS.out");
  ofstream *out[2];
  char *name;
  if(!file) name = "SmoothingSolution_ANALYSIS1.out";
  if( file) name = "SmoothingSolution_ANALYSIS2.out";
  out[file] = new ofstream(name);
  an.Print("\n\n* * * SOLUCAO PELO ANALYSIS: primeira solução  * * *\n\n",*out[file]);
  fBegin = clock();
  mat->SetForcingFunction(0);
  REAL normsol = Norm(Solution());

  while(iter < numiter && normsol < tol) {
    
    fBegin = clock();
    an.Solution().Zero();
    SetDeltaTime(anmesh,mat,iter);
    an.Run();
    CoutTime(fBegin,"TPZNonLinMultGridAnalysis:: Fim system solution actual iteration");
    CoutTime(fInit,"TPZNonLinMultGridAnalysis:: accumulated time");
    fBegin = clock();
    an.LoadSolution();
    cout << "TPZNonLinMultGridAnalysis::SmoothingSolution iteracao = " << ++iter << endl;
    normsol = Norm(Solution());
  }
  if(iter < numiter){
    cout << "\nTPZNonLinMultGridAnalysis::SmoothingSolution the iterative process stopped"
	 << " due the great norm of the solution, norm solution = " << normsol << endl;
  }
  an.LoadSolution();
  an.Print("\n\n* * * SOLUCAO PELO ANALYSIS: última solução  * * *\n\n",*out[file]);
  CoutTime(fInit,"TPZNonLinMultGridAnalysis::SmoothingSolution general time of iterative process");
  out[file]->flush();
  out[file]->close();
}

void TPZNonLinMultGridAnalysis::SmoothingSolution2(REAL tol,int numiter,TPZMaterial *mat,TPZAnalysis &an,
						   ofstream &dxout,int marcha) {

  TPZVec<char *> scalar(1),vector(0);
  scalar[0] = "pressure";
  int dim = mat->Dimension();
  TPZCompMesh *anmesh = an.Mesh();
  ResetReference(anmesh);//retira referências para criar graph consistente
  TPZDXGraphMesh graph(anmesh,dim,mat,scalar,vector);
  SetReference(anmesh);//recupera as referências retiradas
  //ofstream dxout = new ofstream("EulerConsLaw.dx");
  cout << "\nTPZNonLinMultGridAnalysis::IterativeProcess out file : EulerConsLaw.dx\n";
  graph.SetOutFile(dxout);
  int resolution = 0;
  graph.SetResolution(resolution);
  graph.DrawMesh(dim);
  int iter = 0,draw=0;
  fInit = clock();
  an.Solution().Zero();
  fBegin = clock();
  SetDeltaTime(anmesh,mat,iter);
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
    SetDeltaTime(anmesh,mat,iter);
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
      dxout.flush();
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
  dxout.flush();
  ofstream out("SmoothingSolution2_ANALYSIS.out");
  an.Print("\n\n* * * SOLUCAO PELO ANALYSIS: última solução  * * *\n\n",out);
  out.flush();
  out.close();
}

void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix &sol,TPZAnalysis &an,
					     char *decompose,TPZFMatrix &res){

  //TPZStepSolver *solver = dynamic_cast<TPZStepSolver *>(&an.Solver());
  //TPZMatrix *stiff = solver->Matrix();
  TPZMatrix *stiff = an.Solver().Matrix();
  ofstream out("CalcResidual_STIFF.out");
  stiff->Print("\n\n\t\t\t* * * MATRIZ DE RIGIDEZ * * *\n\n",out);
  int dim = stiff->Dim(),i,j;
  res.Redim(dim,1);
  //cálculo de stiff * solution
  TPZFMatrix tsup(dim,1),diag(dim,1),tinf(dim,1);

  if( !strcmp(decompose , "LDLt") ){
    //triângulo superior
    for(i=0;i<dim;i++){
      REAL sum = 0.;
      for(j=i+1;j<dim;j++){
	 sum += stiff->GetVal(i,j) * sol(j,0);
      }
      tsup(i,0) = sol(i,0) + sum;
    }
    //diagonal
    for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,0);
    //triângulo superior
    for(i=0;i<dim;i++){
      REAL sum = 0.;
      for(j=0;j<i;j++){
	 sum += stiff->GetVal(i,j) * diag(j,0);
      }
      tinf(i,0) = sum + sol(i,0);
    }
    //diferenca (f - stiff * x)
    TPZFMatrix rhs = an.Rhs();
    for(i=0;i<dim;i++) res(i,0) = rhs.GetVal(i,0) - tinf(i,0);
    ofstream out("MATRIZES_DECOMPOSIÇÂO.out");
    sol.Print("\n* * * sol * * *\n",out);
    tsup.Print("\n* * * tsup * * *\n",out);
    diag.Print("\n* * * diag * * *\n",out);
    tinf.Print("\n* * * tinf * * *\n",out);
    rhs.Print("\n* * * rhs * * *\n",out);
    res.Print("\n* * * res * * *\n",out);
    stiff->Print("\n* * * stiff * * *\n",out);
    cout << "\nNorma do resíduo: " << Norm(res) << endl;
  } else {
    cout << "TPZNonLinMultGridAnalysis::CalcResidual Calculation of the residue for this"
	 << " decomposition is not implemented, implements now!\n";
  }
}


void TPZNonLinMultGridAnalysis::CalcResidual(TPZMatrix &sol,TPZFMatrix &anres,TPZFMatrix &res,
					     TPZAnalysis &an,char *decompose){

  //TPZStepSolver *solver = dynamic_cast<TPZStepSolver *>(&an.Solver());
  //TPZMatrix *stiff = solver->Matrix();
  TPZMatrix *stiff = an.Solver().Matrix();
  ofstream out("CalcResidual_STIFF.out");
  stiff->Print("\n\n\t\t\t* * * MATRIZ DE RIGIDEZ * * *\n\n",out);
  int dim = stiff->Dim(),i,j;
  //cálculo de stiff * solution
  TPZFMatrix tsup(dim,1),diag(dim,1),tinf(dim,1);

  if( !strcmp(decompose , "LDLt") ){
    //triângulo superior
    for(i=0;i<dim;i++){
      REAL sum = 0.;
      for(j=i+1;j<dim;j++){
	 sum += stiff->GetVal(i,j);
      }
      tsup(i,0) = sol(i,0) + sum;
    }
    //diagonal
    for(i=0;i<dim;i++) diag(i,0) = stiff->GetVal(i,i) * tsup(i,0);
    //triângulo superior
    for(i=0;i<dim;i++){
      REAL sum = 0.;
      for(j=0;j<i-1;j++){
	 sum += stiff->GetVal(i,j) * diag(i,0);
      }
      tinf(i,0) = sum + sol(i,0);
    }
    //diferenca (f - stiff * x)
    //TPZFMatrix rhs = an.Rhs();
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

void TPZNonLinMultGridAnalysis::TwoGridAlgorithm(ostream &out,int nummat){

  TPZCompMesh *coarcmesh = fMeshes[0];//malha grosseira inicial
  int meshdim = coarcmesh->Dimension();
  //criando a malha fina
  int levelnumbertorefine = 1;
  cout << "TPZNonLinMultGridAnalysis:: número de níveis a dividir: ";
  //cin >> levelnumbertorefine;
  int setdegree = -1;//preserva o grau da malha inicial
  //newmesh = 0: coarcmesh se tornou a malha fina
  TPZCompMesh *finemesh = UniformlyRefineMesh(coarcmesh,levelnumbertorefine,setdegree);
  finemesh->SetDimModel(meshdim);
  finemesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL FINA * * *\n\n");
  //obtendo-se a malha menos fina por agrupamento
  int levelnumbertogroup = 0;//sera obtido por agrupamento o nível 0
  cout << "TPZNonLinMultGridAnalysis:: número do nível a ser agrupado: ";
  //cin >> levelnumbertogroup;
  TPZCompMesh *aggmesh = AgglomerateMesh(finemesh,levelnumbertogroup);
  aggmesh->SetName("\n\t\t\t* * * MALHA COMPUTACIONAL AGLOMERADA * * *\n\n");
  aggmesh->SetDimModel(meshdim);
  AppendMesh(aggmesh);
  aggmesh->Reference()->SetName("\n\t\t\t* * * MALHA GEOMÈTRICA REFINADA * * *\n\n");
  aggmesh->Reference()->Print(out);//malha geométrica é uma só
  //out << "\n\n\t\t\t* * * MALHA AGLOMERADA * * *\n\n";
  aggmesh->Print(out);
  //out << "\n\n\t\t\t* * * MALHA FINA * * *\n\n";
  finemesh->Print(out);
  out.flush();
  //analysis na malha aglomerada
  TPZAnalysis coarsean(fMeshes[1]);
  TPZSkylineStructMatrix coarsestiff(fMeshes[1]);
  coarsean.SetStructuralMatrix(coarsestiff);
  TPZStepSolver coarsesolver;
  coarsesolver.SetDirect(ELDLt);
  TPZMatrixSolver *clone = dynamic_cast<TPZMatrixSolver *>(coarsesolver.Clone());
  coarsean.SetSolver(*clone);
  fSolvers.Push(clone);
  coarsean.Solution().Zero();
  fSolutions.Push(new TPZFMatrix(coarsean.Solution()));
  fPrecondition.Push(0);
  //analysis na malha fina
  AppendMesh(finemesh);
  TPZAnalysis finean(fMeshes[2]);
  TPZSkylineStructMatrix finestiff(fMeshes[2]);
  finean.SetStructuralMatrix(finestiff);
  TPZStepSolver finesolver;
  finesolver.SetDirect(ELDLt);
  clone = dynamic_cast<TPZMatrixSolver *>(finesolver.Clone());
  finean.SetSolver(*clone);
  fSolvers.Push(clone);
  finean.Solution().Zero();
  fSolutions.Push(new TPZFMatrix(finean.Solution()));
  fPrecondition.Push(0);
  //suavisar a solução na malha fina
  REAL tol = 1.e15;//valor máximo da ||solução||
  int numiter,marcha;
  cout << "\nNumero de iteracoes pre-suavisamento? :\n";
  cin >> numiter;
  //numiter = 10;
  cout << "main:: Parametro marcha (nulo não cria .dx): 0\n";
  //cin >> marcha;
  marcha = 0;
  fMeshes[2]->Reference()->ResetReference();
  fMeshes[2]->LoadReferences();
  TPZMaterial *mat = fMeshes[2]->FindMaterial(nummat);
  ofstream *dxout = new ofstream("PreSmoothingEuler.dx");
  SmoothingSolution(tol,numiter,mat,finean,*dxout,marcha);
  finean.Print("\n\n\t\t\t* * * ANALYSIS MALHA FINA após SmoothingSolution * * *\n\n",out);
  out.flush();
  // TRANSFERÊNCIA DE SOLUÇÕES
  ofstream OUT("TwoGridAlgorithm_TRANSFER.out");
  TPZTransfer transfer;
  fMeshes[2]->BuildTransferMatrixDesc(*fMeshes[1],transfer);
  TPZFMatrix projectsol;
  //fMeshes[1]->ProjectSolution(projectsol);
  //return;
  int neq;// = fMeshes[1]->NEquations();

  if(1){//TESTE 3
    cout << "\n\nTESTE 3: Project Solution\n\n";
    //ESTE TESTE MOSTRA QUE A SOLUÇÃO É TRANSFERIDA DA 
    //MALHA FINA PARA A MALHA GROSSA (EXATA PARA GRAU 0)
    neq = fMeshes[2]->NEquations();
    TPZFMatrix finesol(neq,1,1.);
    REAL val0 = 1.,val1 = 0.,val2 = 0.1,val3 = 1.0;//com CFL = 0
    //cout << "TPZNonLinMultGridAnalysis::TwoGridAlgorithm: Entre valores : ";
    //cin >> val1 >> val2 >> val3;
    REAL eps = 0.0;
    for(int i=0;i<neq;i++){
      int k = (i+1)%4;
      if(k == 2){
	finesol(i,0) = val1;// + eps * i;
	continue;
      }
      if(k == 3){
	finesol(i,0) = val2 + eps * i;
	continue;
      }
      if(k == 0) finesol(i,0) = val3 + eps * i;
      if(k == 1) finesol(i,0) = val0 + eps * i;
    }
    fMeshes[2]->LoadSolution(finesol);
    neq = fMeshes[1]->NEquations();
    TPZFMatrix coarsesol(neq,1,0.);
    fMeshes[1]->ProjectSolution(coarsesol);
//     TPZFMatrix DIFF = fMeshes[1]->Solution() - coarsesol;
//     cout << "* * * NORMA DO ERROR -->>>" << Norm(DIFF) << endl;
    fMeshes[1]->LoadSolution(coarsesol);
    //transfer.TransferSolution(coarsesol,finesol);
    finesol.Print("\n\n\t\t\t* * * finesol * * *\n\n",OUT);
    coarsesol.Print("\n\n\t\t\t* * * coarsesol * * *\n\n",OUT);
    int nel = fMeshes[1]->ElementVec().NElements(),i;
    TPZVec<REAL> csi(3,0.),x(3);
    TPZAgglomerateElement *aggel;
    for(i=0;i<nel;i++){
      TPZCompEl *cel = fMeshes[1]->ElementVec()[i];
      if(!cel) continue;
      if(cel->Type() != EAgglomerate) continue;
      aggel = dynamic_cast<TPZAgglomerateElement *>(cel);
//       TPZStack<int> elvec;
//       aggel->IndexesDiscSubEls(elvec);
//       TPZCompEl *fineel0 = fMeshes[2]->ElementVec()[elvec[0]];//basta 1 para teste
      int nindex = aggel->NIndexes(),i;
      for(i=0;i<nindex;i++){
	TPZCompEl *fineeli = aggel->FineElement(i);
	TPZGeoEl *ref = fineeli->Reference();
	TPZManVector<REAL> sol;
	int var = 4;//pressão
	fineeli->Solution(csi,var,sol);
	ref->X(csi,x);
	OUT << "Elemento -> " << ref->Id() << "\n";
	OUT << "x = " << x[0] << " " << x[1] << " " << x[2] << " : ";
	OUT << "solução fina -> " << sol[0] << "\n";
	aggel->Solution(x,-var,sol);//ponto real
	OUT << "\t\t\tsolução grossa -> " << sol[0] << "\n";
      }
    }
    OUT.flush();
    OUT.close();
  }

  if(0){//TESTE 2
  cout << "\n\nTESTE 2: Transfer Solution\n\n";
  //ESTE TESTE MOSTRA QUE A SOLUÇÃO É TRANSFERIDA DA 
  //MALHA GROSSA PARA A MALHA FINA DE FORMA EXATA
  //PARA QUALQUER GRAU DE INTERPOLAÇÃO
  neq = fMeshes[1]->NEquations();
  TPZFMatrix coarsesol(neq,1,1.);
  REAL val0 = 1.,val1 = 0.,val2 = 0.1,val3 = 1.0;//com CFL = 0
  //cout << "TPZNonLinMultGridAnalysis::TwoGridAlgorithm: Entre valores : ";
  //cin >> val1 >> val2 >> val3;
  REAL eps = 0.001;
  for(int i=0;i<neq;i++){
    int k = (i+1)%4;
    if(k == 2){
      coarsesol(i,0) = val1;// + eps * i;
      continue;
    }
    if(k == 3){
      coarsesol(i,0) = val2 + eps * i;
      continue;
    }
    if(k == 0) coarsesol(i,0) = val3 + eps * i;
    if(k == 1) coarsesol(i,0) = val0 + eps * i;
  }
  neq = fMeshes[2]->NEquations();
  TPZFMatrix finesol(neq,1,0.);
  fMeshes[1]->LoadSolution(coarsesol);
  transfer.TransferSolution(coarsesol,finesol);
  coarsesol.Print("\n\n\t\t\t* * * coarsesol * * *\n\n",OUT);
  finesol.Print("\n\n\t\t\t* * * finesol * * *\n\n",OUT);
  fMeshes[2]->LoadSolution(finesol);
  int nel = fMeshes[1]->ElementVec().NElements(),i;
  TPZVec<REAL> csi(3,0.),x(3);
  TPZAgglomerateElement *aggel;
  for(i=0;i<nel;i++){
    TPZCompEl *cel = fMeshes[1]->ElementVec()[i];
    if(!cel) continue;
    if(cel->Type() != EAgglomerate) continue;
    aggel = dynamic_cast<TPZAgglomerateElement *>(cel);
//     TPZStack<int> elvec;
//     aggel->IndexesDiscSubEls(elvec);
//     //basta 1 para teste
//     TPZCompEl *fineel0 = fMeshes[2]->ElementVec()[elvec[0]];
    TPZCompEl *fineel0 = aggel->FineElement(0);//basta 1 para teste
    TPZGeoEl *ref = fineel0->Reference();
    TPZManVector<REAL> sol;
    int var = 4;//pressão
    fineel0->Solution(csi,var,sol);
    ref->X(csi,x);
    OUT << "Elemento -> " << ref->Id() << "\n";
    OUT << "x = " << x[0] << " " << x[1] << " " << x[2] << " : ";
    OUT << "solução fina -> " << sol[0] << "\n";
    aggel->Solution(x,-var,sol);//ponto real
    OUT << "\t\t\tsolução grossa -> " << sol[0] << "\n";
  }
  OUT.flush();
  OUT.close();
  }

  if(0){//TESTE 1: para verificar a qualidade da 
        //         transferência do resíduo
  neq = fMeshes[1]->NEquations();
  TPZFMatrix coarsesol(neq,1,1.),transcoarsesol(neq,1,0.);
  neq = fMeshes[2]->NEquations();
  TPZFMatrix finesol(neq,1,0.);
  transfer.TransferSolution(coarsesol,finesol);
  transfer.TransferResidual(finesol,transcoarsesol);
  coarsesol.Print("\n\n\t\t\t* * * coarsesol * * *\n\n",OUT);
  finesol.Print("\n\n\t\t\t* * * finesol * * *\n\n",OUT);
  transcoarsesol.Print("\n\n\t\t\t* * * transcoarsesol * * *\n\n",OUT);
  OUT.flush();
  OUT.close();
  }

  if(0){
  transfer.Print("\n\n\t\t\t* * * MATRIZ DE TRANSFERÊNCIA * * *\n\n",OUT);
  int neq = fMeshes[1]->NEquations();
  TPZFMatrix coarseRhs(neq,1,0.);
  transfer.TransferResidual(finean.Rhs(),coarseRhs);//transferindo résiduo
  coarseRhs.Print("\n\n* * * RESÍDUO TRANSFERIDO: FINO -> GROSSO * * *\n\n",OUT);
  ofstream *dxout2 = new ofstream("PostSmoothingEuler.dx");
  cout << "\nNumero de iteracoes pós-suavisamento? :\n";
  cin >> numiter;
  fMeshes[1]->Reference()->ResetReference();
  fMeshes[1]->LoadReferences();
  mat = fMeshes[1]->FindMaterial(nummat);
  SmoothingSolution(tol,numiter,mat,coarsean,*dxout2,marcha);
  coarsean.Rhs().Print("\n\n* * * RESÍDUO DA MALHA GROSSA * * *\n\n",OUT);
  coarsean.Print("\n\n\t\t\t* * * ANALYSIS MALHA GROSSA após SmoothingSolution * * *\n\n",out);
  neq = fMeshes[2]->NEquations();
  TPZFMatrix finesol(neq,1,0.);
  transfer.TransferSolution(coarsean.Solution(),finesol);//transferindo solução
  fMeshes[2]->LoadSolution(finesol);
  TPZFMatrix res;
  CalcResidual(finesol,finean,"LDLt",res);
  res.Print("\n\n\t\t\t* * * RESÍDUO MALHA FINA * * *\n\n",out);
  finesol.Print("\n\n* * * SOLUÇÂO TRANSFERIDA: GROSSA -> FINA * * *\n\n",OUT);
  OUT.flush();
  OUT.close();
  }
}

//  CalcResidual(finean.Solution(),res,finean,"LDLt");
//  res.Print("\n\n\t\t\t* * * RESÍDUO MALHA FINA * * *\n\n",out);
