// -*- c++ -*-
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzconnect.h"
#include "pzgmesh.h"
#include "pzfmatrix.h"
#include "pzcompel.h"
#include "TPZConservationLaw.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZCompElDisc.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzelmat.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzmanvector.h"
#include "pzv3dmesh.h"
#include "pzdxmesh.h"
#include "pzmvmesh.h"

#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzmetis.h"
//#include "pzsloan.h"

#include "pzmaterial.h"

#include "pzstrmatrix.h"

#include <fstream>
using namespace std;
#include <string.h>
#include <stdio.h>
using namespace std;



void TPZAnalysis::SetStructuralMatrix(TPZStructMatrix &strmatrix){
    if(fStructMatrix) delete fStructMatrix;
    fStructMatrix = strmatrix.Clone();
}
TPZAnalysis::TPZAnalysis() : fRhs(), fSolution(), fTable() {
	fCompMesh = 0;
	fGeoMesh = 0;
	fStructMatrix = 0;
   fGraphMesh[0] = 0;
   fGraphMesh[1] = 0;
   fGraphMesh[2] = 0;
   fSolver = 0;
}


TPZAnalysis::TPZAnalysis(TPZCompMesh *mesh,std::ostream &out) :
    fRhs(), fSolution(), fTable()
{

	fCompMesh = mesh;
   fGeoMesh = mesh->Reference();
   fStructMatrix = 0;
   fGraphMesh[0] = 0;
   fGraphMesh[1] = 0;
   fGraphMesh[2] = 0;
   fSolver = 0;
   SetBlockNumber();
   fStep = 0;
   fTime = 0.;
   fCompMesh->ExpandSolution();
   fSolution = fCompMesh->Solution();
   int neq = fCompMesh->NEquations();
   fSolution.Resize(neq,1);
}

TPZAnalysis::~TPZAnalysis(void){
	if(fStructMatrix) delete fStructMatrix;
	if(fSolver) delete fSolver;
	int dim;
	for(dim=0; dim<3; dim++) {
		if(fGraphMesh[dim]) delete fGraphMesh[dim];
		fGraphMesh[dim] = 0;
	}
}

void TPZAnalysis::SetBlockNumber(){

	if(!fCompMesh) return;
	fCompMesh->InitializeBlock();

	TPZStack<int> elgraph,elgraphindex;
	int nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int nel = elgraphindex.NElements()-1;
	TPZMetis metis(nel,nindep);
	metis.SetElementGraph(elgraph,elgraphindex);
	TPZVec<int> perm,iperm;
//  ofstream write("Graph.chk");
//  metis.Print(out,"Graph File");
//  metis.Print(write);
//  out.flush();
	metis.Resequence(perm,iperm);
//	fCompMesh->Permute(perm);
}

void TPZAnalysis::Assemble() {
	if(!fCompMesh || !fStructMatrix || !fSolver) return;

	fRhs.Redim(fCompMesh->NEquations(),1);
	fSolver->SetMatrix(0);
	fSolver->SetMatrix(fStructMatrix->CreateAssemble(fRhs));//aqui TPZFMatrix não é nula
//fRhs->Print("Rhs");
//cout.flush();
}

void TPZAnalysis::Solve() {
	int numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return;

	TPZFMatrix residual(fRhs);
	TPZFMatrix delu(numeq,1);
	if(fSolution.Rows() != numeq) {
	  fSolution.Redim(numeq,1);
	} else {
	  fSolver->Matrix()->Residual(fSolution,fRhs,residual);
	}
	//      REAL normres  = Norm(residual);
	//	cout << "TPZAnalysis::Solve residual : " << normres << " neq " << numeq << endl;
	fSolver->Solve(residual, delu);
	fSolution += delu;
	
	fCompMesh->LoadSolution(fSolution);
}

void TPZAnalysis::LoadSolution(){

  if(fCompMesh) {
    fCompMesh->LoadSolution(fSolution);
  }
}

void TPZAnalysis::Print( char *name , std::ostream &out )
{

  out<<endl<<name<<endl;
  int i,nelements = fCompMesh->ConnectVec().NElements();
  for(i=0;i<nelements;i++) {
    TPZConnect &gnod = fCompMesh->ConnectVec()[i];
    if(gnod.SequenceNumber()!=-1) {
    	out << "Connect index " << i << endl;
    	gnod.Print(*fCompMesh,out);
    }
  }

	fSolution.Print("fSolution",out);
   fCompMesh->ConnectSolution(out);

}


void
TPZAnalysis::PostProcess(TPZVec<REAL> &, std::ostream &out )
{

  int neq = fCompMesh->NEquations();
  TPZVec<REAL> ux((int) neq);
	TPZVec<REAL> sigx((int) neq);
	TPZVec<REAL> values(3,0.);
	fCompMesh->LoadSolution(fSolution);
//	SetExact(&Exact);
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	double true_error=0.
			,L2_error=0.
			,estimate=0.;
	int nel = elvec.NElements();
	for(int i=0;i<nel;i++) {
		TPZCompEl *el = (TPZCompEl *) elvec[i];
      if(el) {
           el->EvaluateError(fExact,true_error,L2_error,0,estimate);
           values[0]+=true_error;
           values[1]+=L2_error;
           values[2]+=estimate;
      }
	}


		out<<endl<<"############"<<endl;
		out<<endl<<"true_error="<<values[0]<<endl;
		out<<endl<<"L2_error="<<values[1]<<endl;
		out<<endl<<"estimate="<<values[2]<<endl;

	return;
}

void TPZAnalysis::PostProcessTable( TPZFMatrix &,ostream & )//pos,out
{
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	int nel = elvec.NElements();
	for(int i=0;i<nel;i++) {
		//TPZCompEl *el = (TPZCompEl *) elvec[i];
  	   //if(el)	el->PrintTable(fExact,pos,out);
	}
	return;
}

void TPZAnalysis::ShowShape( TPZVec<char *> &scalnames, TPZVec<char *> &vecnames, char *plotfile, std::ostream &) {//1o :TPZConnect  *nod,

   TPZV3DGraphMesh gg(fCompMesh,2,fCompMesh->MaterialVec()[0]);

	ofstream plot(plotfile);
	gg.SetOutFile(plot);
	gg.SetResolution(0);
	gg.DrawMesh(1);
	gg.DrawSolution(0,0.,scalnames,vecnames);

}

void TPZAnalysis::LoadShape(double ,double , int ,TPZConnect* start){
//void TPZAnalysis::LoadShape(double dx,double dy, int numelem,TPZConnect* start){
	Assemble();
   fRhs.Zero();
   fSolution.Zero();

}

void TPZAnalysis::Run(ostream &out){

	Assemble();
	Solve();
}

void TPZAnalysis::DefineGraphMesh(int dim, TPZVec<char *> &scalnames, TPZVec<char *> &vecnames, char *plotfile) {

  int dim1 = dim-1;
  TPZMaterial *mat;
  int nmat,imat;
  nmat = fCompMesh->MaterialVec().NElements();
  for(imat=0; imat<nmat; imat++) {
    mat = fCompMesh->MaterialVec()[imat];
    if(mat && mat->Id() >= 0 && mat->Dimension() == dim) break;
    mat = 0;
  }
  if(fGraphMesh[dim1]) delete fGraphMesh[dim1];
  fScalarNames[dim1] = scalnames;
  fVectorNames[dim1] = vecnames;
  //misael
  if(! ( strcmp( strrchr(plotfile,'.')+1,"plt") ) )	{
    fGraphMesh[dim1] = new TPZV3DGraphMesh(fCompMesh,dim,mat) ;
  }else if(!strcmp( strrchr(plotfile,'.')+1,"dx") ) {
    fGraphMesh[dim1] = new TPZDXGraphMesh(fCompMesh,dim,mat,scalnames,vecnames) ;
  }else if(!strcmp( strrchr(plotfile,'.')+1,"pos") ) {
    fGraphMesh[dim1] = new TPZMVGraphMesh(fCompMesh,dim,mat);
  } else {
    cout << "grafgrid was not created\n";
    fGraphMesh[dim1] = 0;
  }
  if(fGraphMesh[dim1]) {
    ofstream *plot = new ofstream(plotfile);
    fGraphMesh[dim1]->SetOutFile(*plot);
  }
}

void TPZAnalysis::PostProcess(int resolution) {
	int dim;
   for(dim=1; dim<=3; dim++) {
   	PostProcess(resolution, dim);
   }
}

void TPZAnalysis::PostProcess(int resolution, int dimension){
  int dim1 = dimension-1;
  if(!fGraphMesh[dim1]) return;
  TPZMaterial *mat;
  int nmat,imat;
  nmat = fCompMesh->MaterialVec().NElements();
  for(imat=0; imat<nmat; imat++) {
    mat = fCompMesh->MaterialVec()[imat];
    if(mat && mat->Id() >= 0 && mat->Dimension() == dimension) break;
    mat = 0;
  }
  if(!mat) return;
  fGraphMesh[dim1]->SetCompMesh(fCompMesh,mat);

  fGraphMesh[dim1]->SetResolution(resolution);
  fGraphMesh[dim1]->DrawMesh(1);
  fGraphMesh[dim1]->DrawSolution(fStep,fTime);
//   delete fGraphMesh[dim1];
//   fGraphMesh[dim1] = 0;
}

void TPZAnalysis::AnimateRun(int num_iter, int steps,
			     TPZVec<char *> &scalnames, TPZVec<char *> &vecnames, char *plotfile)
{
  Assemble();
  
  int numeq = fCompMesh->NEquations();
  if(fRhs.Rows() != numeq ) return;
  
  TPZFMatrix residual(fRhs);
  int dim = HighestDimension();
  TPZMaterial *mat;
  int nmat,imat;
  nmat = fCompMesh->MaterialVec().NElements();
  for(imat=0; imat<nmat; imat++) {
    mat = fCompMesh->MaterialVec()[imat];
    if(mat && mat->Id() >= 0 && mat->Dimension() == dim) break;
    mat = 0;
  }
  
  TPZDXGraphMesh gg(fCompMesh,dim,mat,scalnames,vecnames) ;
  // 		TPZV3DGrafGrid gg(fCompMesh) ;
  ofstream plot(plotfile);
  gg.SetOutFile(plot);
  gg.SetResolution(0);
  gg.DrawMesh(num_iter);
  
  int i;
  for(i=1; i<=num_iter;i+=steps){
    
 
    TPZStepSolver sol;
    sol.ShareMatrix(Solver());
    sol.SetJacobi(i,0.,0);
    SetSolver(sol);
    //    Solver().SetNumIterations(i);
    fSolver->Solve(fRhs, fSolution);
    
    fCompMesh->LoadSolution(fSolution);
    gg.DrawSolution(i-1,0);
  }
}

int TPZAnalysis::HighestDimension(){
	TPZAdmChunkVector<TPZMaterial *> &matmap = fCompMesh->MaterialVec();
   int nel = matmap.NElements();
   int i,dim = 0;
   for(i=0;i<nel;i++) {
   	TPZMaterial *mat = (TPZMaterial *) matmap[i];
      if(mat) dim = mat->Dimension() > dim ? mat->Dimension() : dim;
   }
   return dim;
}

TPZAnalysis::TTablePostProcess::TTablePostProcess() :
	fGeoElId(), fCompElPtr(), fLocations(), fVariableNames() {
   fDimension = -1;
   fOutfile = 0;
}

TPZAnalysis::TTablePostProcess::~TTablePostProcess() {
	fDimension = -1;
   int numvar = fVariableNames.NElements();
   int iv;
   for(iv=0; iv<numvar; iv++) {
     	char *name = (char *) fVariableNames[iv];
     	if(name) delete name;
   }
   if(fOutfile) delete fOutfile;
   fOutfile = 0;
}

void TPZAnalysis::DefineElementTable(int dimension, TPZVec<int> &GeoElIds, TPZVec<REAL> &points) {
	fTable.fDimension = dimension;
   fTable.fGeoElId = GeoElIds;
   fTable.fLocations = points;
   int numel = GeoElIds.NElements();
   fTable.fCompElPtr.Resize(numel);
   int iel;
   for(iel=0; iel<numel; iel++) {
   	TPZGeoEl *gel = (TPZGeoEl *) fGeoMesh->FindElement(GeoElIds[iel]);
   	fTable.fCompElPtr[iel] = (gel) ? gel->Reference() : 0;
   }
}

void TPZAnalysis::SetTablePostProcessFile(char *filename) {
	if(fTable.fOutfile) delete fTable.fOutfile;
   fTable.fOutfile = new ofstream(filename);
}

void TPZAnalysis::SetTableVariableNames(int numvar, char **varnames) {
   int nvar = fTable.fVariableNames.NElements();
   int iv;
   for(iv=0; iv<nvar; iv++) {
     	char *name = (char *) fTable.fVariableNames[iv];
     	if(name) delete[] name;
   }
   fTable.fVariableNames.Resize(numvar);
   for(iv=0; iv<numvar; iv++) {
   	char *name = new char[strlen(varnames[iv]+1)];
      strcpy(name,varnames[iv]);
      fTable.fVariableNames[iv] = name;
   }
}


void TPZAnalysis::PrePostProcessTable(){
	TPZCompEl *cel;
   int numvar = fTable.fVariableNames.NElements();
   for(int iv=0; iv<numvar; iv++) {
      int numel = fTable.fCompElPtr.NElements();
      for(int iel=0; iel<numel; iel++) {
         cel = (TPZCompEl *) fTable.fCompElPtr[iel];
         if(cel) cel->PrintTitle((char *)fTable.fVariableNames[iv],*(fTable.fOutfile));
      }
   }
   *(fTable.fOutfile) << endl;
   int dim;
   TPZVec<REAL> point(fTable.fDimension);
   for(dim=1; dim<fTable.fDimension+1; dim++) {
      for(int iv=0; iv<numvar; iv++) {
         int numel = fTable.fCompElPtr.NElements();
         for(int iel=0; iel<numel; iel++) {
         	int d;
            for(d=0; d<fTable.fDimension; d++) {
            	point[d] = fTable.fLocations[iel*fTable.fDimension+d];
            }
            cel = (TPZCompEl *) fTable.fCompElPtr[iel];
            if(cel) cel->PrintCoordinate(point,dim,*(fTable.fOutfile));
         }
      }
   *(fTable.fOutfile) << endl;
   }
}

void TPZAnalysis::PostProcessTable() {
   TPZVec<REAL> point(fTable.fDimension);
   int numvar = fTable.fVariableNames.NElements();
   TPZCompEl *cel;
   for(int iv=0; iv<numvar; iv++) {
      int numel = fTable.fCompElPtr.NElements();
      for(int iel=0; iel<numel; iel++) {
         int d;
         for(d=0; d<fTable.fDimension; d++) {
            point[d] = fTable.fLocations[iel*fTable.fDimension+d];
         }
         cel = (TPZCompEl *) fTable.fCompElPtr[iel];
         if(cel) cel->PrintSolution(point,(char*)fTable.fVariableNames[iv],*(fTable.fOutfile));
      }
   }
   *(fTable.fOutfile) << endl;
}
void TPZAnalysis::SetSolver(TPZMatrixSolver &solver){
  if(fSolver) delete fSolver;
    fSolver = (TPZMatrixSolver *) solver.Clone();
}

void CoutTime(clock_t &start,char *title);
void SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat);
static clock_t begin,init;
void TPZAnalysis::IterativeProcess(ostream &out,REAL tol,int numiter,TPZMaterial *mat,int marcha,int resolution) {

  cout << "PZAnalysis::IterativeProcessTest beginning of the iterative process, general time 0\n";
  TPZVec<char *> scalar(1),vector(0);
  scalar[0] = "pressure";
  //scalar[1] = "density";
  //scalar[2] = "normvelocity";
  cout << "TPZAnalysis::IterativeProcess solution required : " << scalar[0] << endl;
  //       << "\n" << scalar[1] << "\n" << scalar[2] << endl;
  int dim = mat->Dimension();
  TPZDXGraphMesh graph(Mesh(),dim,mat,scalar,vector);
  ofstream *dxout = new ofstream("ConsLaw.dx");
  cout << "\nTPZAnalysis::IterativeProcess out file : ConsLaw.dx\n";
  graph.SetOutFile(*dxout);
  graph.SetResolution(resolution);
  graph.DrawMesh(dim);
  int iter = 0,draw=0;
  fSolution.Zero();
  begin = clock();
  Run();
  cout << "TPZAnalysis::IterativeProcess iteracao = " << ++iter << endl;
  CoutTime(begin,"TPZAnalysis:: Fim system solution first iteration");
  begin = clock();
  LoadSolution();
  //CoutTime(begin,"TPZAnalysis:: Fim Load Solution");
  SetDeltaTime(fCompMesh,mat);
  REAL time = (dynamic_cast<TPZConservationLaw *> (mat))->TimeStep();
  begin = clock();
  graph.DrawSolution(draw++,time);
  dxout->flush();
  //CoutTime(begin,"TPZAnalysis:: Fim Draw Solution");
  mat->SetForcingFunction(0);
  REAL normsol = Norm(fSolution);

  while(iter < numiter && normsol < tol) {
    
    begin = clock();
    fSolution.Zero();
    Run();
    CoutTime(begin,"TPZAnalysis:: Fim system solution actual iteration");
    CoutTime(init,"TPZAnalysis:: accumulated time");
    begin = clock();
    LoadSolution();
    //CoutTime(begin,"TPZAnalysis:: Fim Load Solution");
    SetDeltaTime(fCompMesh,mat);
    time = (dynamic_cast<TPZConservationLaw *> (mat))->TimeStep();
    if( REAL(iter) / REAL(marcha) == draw || marcha == 1){
      begin = clock();
      graph.DrawSolution(draw++,time);
      dxout->flush();
      //CoutTime(begin,"TPZAnalysis:: Fim Draw Solution");
    }
    cout << "TPZAnalysis::IterativeProcess iteracao = " << ++iter << endl;
    normsol = Norm(fSolution);
  }
  out.flush();
  dxout->flush();
  if(iter < numiter){
    cout << "\nTPZAnalysis::IterativeProcess the iterative process stopped due the great norm "
	 << "of the solution, norm solution = " << normsol << endl;
  }
  CoutTime(init,"TPZAnalysis:: general time of iterative process");
}

void TPZAnalysis::IterativeProcessTest(ostream &out,REAL tol,int numiter,TPZMaterial *mat,int marcha,int resolution) {

  cout << "PZAnalysis::IterativeProcessTest beginning of the iterative process, general time 0\n";
  TPZVec<char *> scalar(1),vector(0);
  scalar[0] = "Solution";
  cout << "TPZAnalysis::IterativeProcess solution required : " << scalar[0] << endl;
  int dim = mat->Dimension();
  TPZDXGraphMesh graph(Mesh(),dim,mat,scalar,vector);
  ofstream *dxout = new ofstream("ConsLaw.dx");
  cout << "\nTPZAnalysis::IterativeProcess out file : ConsLaw.dx\n";
  graph.SetOutFile(*dxout);
  graph.SetResolution(resolution);
  graph.DrawMesh(dim);
  int iter = 0,draw=0;
  fSolution.Zero();
  begin = clock();
  Run();
  cout << "TPZAnalysis::IterativeProcess iteracao = " << ++iter << endl;
  CoutTime(begin,"TPZAnalysis:: Fim system solution first iteration");
  begin = clock();
  LoadSolution();
  //CoutTime(begin,"TPZAnalysis:: Fim Load Solution");
  begin = clock();
  REAL time = (dynamic_cast<TPZConservationLaw *> (mat))->TimeStep();
  graph.DrawSolution(draw++,time);
  dxout->flush();
  //CoutTime(begin,"TPZAnalysis:: Fim Draw Solution");
  mat->SetForcingFunction(0);
  REAL normsol = Norm(fSolution);

  while(iter < numiter && normsol < tol) {
    
    begin = clock();
    fSolution.Zero();
    Run();
    CoutTime(begin,"TPZAnalysis:: Fim system solution actual iteration");
    begin = clock();
    LoadSolution();
    //CoutTime(begin,"TPZAnalysis:: Fim Load Solution");
    if( REAL(iter) / REAL(marcha) == draw || marcha == 1){
      begin = clock();
      time = (dynamic_cast<TPZConservationLaw *> (mat))->TimeStep();
      graph.DrawSolution(draw++,time);
      dxout->flush();
      //CoutTime(begin,"TPZAnalysis:: Fim Draw Solution");
    }
    cout << "TPZAnalysis::IterativeProcess iteracao = " << ++iter << endl;
    normsol = Norm(fSolution);
  }
  out.flush();
  dxout->flush();
  if(iter < numiter){
    cout << "\nTPZAnalysis::IterativeProcess the iterative process stopped due the great norm "
	 << "of the solution, norm solution = " << normsol << endl;
  }
  CoutTime(init,"TPZAnalysis:: general time of iterative process");
}

void CoutTime(clock_t &start,char *title){
    clock_t end = clock();
    cout << title <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

void SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat){

  int nstate = mat->NStateVariables();
  REAL maxveloc = CompMesh->MaxVelocityOfMesh(nstate);
  //REAL deltax = CompMesh->DeltaX();
  REAL deltax = CompMesh->LesserEdgeOfMesh();
  //REAL deltax = CompMesh->MaximumRadiusOfEl();
  TPZCompElDisc *disc;
  int degree = disc->gDegree;
  REAL CFL = 1./((2.0*(REAL)degree) + 1.0);
  TPZDiffusionConsLaw *diff;
  //  if(!degree) CFL = diff->fCFL;

  dynamic_cast<TPZConservationLaw *>(mat)->SetTimeStep(CFL*deltax/maxveloc);
}

/*
  int sizesol = fSolution.Rows();

  if(0 && sizesol < 1000 && iter < 20){
    out << "TPZAnalysis::IterativeProcess iteracao = " << iter 
	<< " : norma da solucao ||Solution||: " << Norm(fSolution) << endl;
    Print("FEM SOLUTION ",out);
  }
*/
/*
    if(0 && sizesol < 1000 && iter < 20){
      out << "TPZAnalysis::IterativeProcess iteracao = " << iter 
	  << " : norma da solucao ||Solution||: " << Norm(fSolution) << endl;
      Print("FEM SOLUTION ",out);
    }
*/
