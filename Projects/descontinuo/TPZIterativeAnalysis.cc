// -*- c++ -*-
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "TPZFlowCMesh.h"
#include "TPZConservationLaw.h"
#include "TPZIterativeAnalysis.h"
#include "TPZAgglomerateEl.h"
#include <string.h>
#include <stdio.h>
#include "pzdxmesh.h"
using namespace std;

TPZIterativeAnalysis::TPZIterativeAnalysis(TPZCompMesh *mesh,std::ostream &out) : TPZAnalysis(mesh,out), fBegin(0),fInit(0) {

  //fBegin = 0.0;
  //fInit = 0.0;
}

void TPZIterativeAnalysis::IterativeProcess(ostream &out,REAL tol,int numiter,TPZMaterial *mat,int marcha,int resolution) {

  cout << "PZAnalysis::IterativeProcessTest beginning of the iterative process, general time 0\n";
  TPZVec<char *> scalar(1),vector(0);
  scalar[0] = "pressure";
  //scalar[1] = "density";
  //scalar[2] = "normvelocity";
  cout << "TPZIterativeAnalysis::IterativeProcess solution required : " << scalar[0] << endl;
  //       << "\n" << scalar[1] << "\n" << scalar[2] << endl;
  int dim = mat->Dimension();
  ResetReference(Mesh());//retira referências para criar graph consistente
  TPZDXGraphMesh graph(Mesh(),dim,mat,scalar,vector);
  SetReference(Mesh());//recupera as referências retiradas
  ofstream *dxout = new ofstream("ConsLaw.dx");
  cout << "\nTPZIterativeAnalysis::IterativeProcess out file : ConsLaw.dx\n";
  graph.SetOutFile(*dxout);
  graph.SetResolution(resolution);
  graph.DrawMesh(dim);
  int iter = 0,draw=0;
//   SetReference(Mesh());//recupera as referências retiradas
  fSolution.Zero();
  fBegin = clock();
  Run();
  cout << "TPZIterativeAnalysis::IterativeProcess iteracao = " << ++iter << endl;
  CoutTime(fBegin,"TPZIterativeAnalysis:: Fim system solution first iteration");
  fBegin = clock();
  LoadSolution();
  //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Load Solution");
  SetDeltaTime(fCompMesh,mat);
  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
  REAL time = law->TimeStep();
  fBegin = clock();
  graph.DrawSolution(draw++,time);
  dxout->flush();
  mat->SetForcingFunction(0);
  REAL normsol = Norm(fSolution);

  while(iter < numiter && normsol < tol) {
    
    fBegin = clock();
    fSolution.Zero();
    Run();
    CoutTime(fBegin,"TPZIterativeAnalysis:: Fim system solution actual iteration");
    CoutTime(fInit,"TPZIterativeAnalysis:: accumulated time");
    fBegin = clock();
    LoadSolution();
    //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Load Solution");
    SetDeltaTime(fCompMesh,mat);
    if( REAL(iter) / REAL(marcha) == draw || marcha == 1){
      fBegin = clock();
      time = law->TimeStep();
      graph.DrawSolution(draw++,time);
      dxout->flush();
      //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Draw Solution");
    }
    cout << "TPZIterativeAnalysis::IterativeProcess iteracao = " << ++iter << endl;
    normsol = Norm(fSolution);
  }
  out.flush();
  dxout->flush();
  if(iter < numiter){
    cout << "\nTPZIterativeAnalysis::IterativeProcess the iterative process stopped due the great norm "
	 << "of the solution, norm solution = " << normsol << endl;
    time = law->TimeStep();
    graph.DrawSolution(draw++,time);
    dxout->flush();
  }
  LoadSolution();
  //dxout->close();
  CoutTime(fInit,"TPZIterativeAnalysis:: general time of iterative process");
}

void TPZIterativeAnalysis::IterativeProcessTest(ostream &out,REAL tol,int numiter,TPZMaterial *mat,int marcha,int resolution) {

  cout << "PZAnalysis::IterativeProcessTest beginning of the iterative process, general time 0\n";
  TPZVec<char *> scalar(1),vector(0);
  scalar[0] = "Solution";
  cout << "TPZIterativeAnalysis::IterativeProcess solution required : " << scalar[0] << endl;
  int dim = mat->Dimension();
  TPZDXGraphMesh graph(Mesh(),dim,mat,scalar,vector);
  ofstream *dxout = new ofstream("ConsLaw.dx");
  cout << "\nTPZIterativeAnalysis::IterativeProcess out file : ConsLaw.dx\n";
  graph.SetOutFile(*dxout);
  graph.SetResolution(resolution);
  graph.DrawMesh(dim);
  int iter = 0,draw=0;
  fSolution.Zero();
  fBegin = clock();
  Run();
  cout << "TPZIterativeAnalysis::IterativeProcess iteracao = " << ++iter << endl;
  CoutTime(fBegin,"TPZIterativeAnalysis:: Fim system solution first iteration");
  fBegin = clock();
  LoadSolution();
  //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Load Solution");
  fBegin = clock();
  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
  REAL time = law->TimeStep();
  graph.DrawSolution(draw++,time);
  dxout->flush();
  //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Draw Solution");
  mat->SetForcingFunction(0);
  REAL normsol = Norm(fSolution);

  while(iter < numiter && normsol < tol) {
    
    fBegin = clock();
    fSolution.Zero();
    Run();
    CoutTime(fBegin,"TPZIterativeAnalysis:: Fim system solution actual iteration");
    CoutTime(fInit,"TPZIterativeAnalysis:: accumulated time");
    fBegin = clock();
    LoadSolution();
    //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Load Solution");
    if( REAL(iter) / REAL(marcha) == draw || marcha == 1){
      fBegin = clock();
      time = law->TimeStep();
      graph.DrawSolution(draw++,time);
      dxout->flush();
      //CoutTime(fBegin,"TPZIterativeAnalysis:: Fim Draw Solution");
    }
    cout << "TPZIterativeAnalysis::IterativeProcess iteracao = " << ++iter << endl;
    normsol = Norm(fSolution);
  }
  out.flush();
  dxout->flush();
  if(iter < numiter){
    cout << "\nTPZIterativeAnalysis::IterativeProcess the iterative process stopped due the great norm "
	 << "of the solution, norm solution = " << normsol << endl;
    time = law->TimeStep();
    graph.DrawSolution(draw++,time);
    dxout->flush();
  }
  CoutTime(fInit,"TPZIterativeAnalysis:: general time of iterative process");
}

void TPZIterativeAnalysis::CoutTime(clock_t &start,char *title){
    clock_t end = clock();
    cout << title <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

void TPZIterativeAnalysis::SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat){

  TPZFlowCompMesh *fm  = dynamic_cast<TPZFlowCompMesh *>(CompMesh);//= new TPZFlowCompMesh(CompMesh->Reference());
  REAL maxveloc = fm->MaxVelocityOfMesh(mat->NStateVariables(),1.4);
  REAL deltax = CompMesh->LesserEdgeOfMesh();//REAL deltax = CompMesh->DeltaX();//REAL deltax = CompMesh->MaximumRadiusOfEl();
  TPZCompElDisc *disc;
  int degree = disc->gDegree;
  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
  law->SetDeltaTime(maxveloc,deltax,degree);
}

/*
  int sizesol = fSolution.Rows();

  if(0 && sizesol < 1000 && iter < 20){
    out << "TPZIterativeAnalysis::IterativeProcess iteracao = " << iter 
	<< " : norma da solucao ||Solution||: " << Norm(fSolution) << endl;
    Print("FEM SOLUTION ",out);
  }
*/
/*
    if(0 && sizesol < 1000 && iter < 20){
      out << "TPZIterativeAnalysis::IterativeProcess iteracao = " << iter 
	  << " : norma da solucao ||Solution||: " << Norm(fSolution) << endl;
      Print("FEM SOLUTION ",out);
    }
*/

void TPZIterativeAnalysis::ResetReference(TPZCompMesh *aggcmesh){
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

void TPZIterativeAnalysis::SetReference(TPZCompMesh *aggcmesh){

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
