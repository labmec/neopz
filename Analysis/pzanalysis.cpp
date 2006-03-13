//$Id: pzanalysis.cpp,v 1.24 2006-03-13 11:54:00 phil Exp $

// -*- c++ -*-
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzconnect.h"
#include "pzgmesh.h"
#include "pzfmatrix.h"
#include "pzcompel.h"
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
#include "pzsloan.h"
#include "pzmaterial.h"
#include "pzstrmatrix.h"

#include "tpznodesetcompute.h"
#include "tpzsparseblockdiagonal.h"
#include "pzseqsolver.h"
#include "pzbdstrmatrix.h"

#include "pzlog.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.analysis"));
#endif

#include <fstream>
using namespace std;
#include <string>
#include <stdio.h>
#include <stdlib.h>
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
	this->SetCompMesh(mesh);
}

void TPZAnalysis::SetCompMesh(TPZCompMesh * mesh) {
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
	TPZVec<int> perm,iperm;

	TPZStack<int> elgraph,elgraphindex;
	int nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int nel = elgraphindex.NElements()-1;
	//TPZSloan sloan(nel,nindep);
	//sloan.SetElementGraph(elgraph,elgraphindex);
//	TPZVec<int> perm,iperm;
	//sloan.Resequence(perm,iperm);
	//fCompMesh->Permute(perm);
/*
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);

	TPZMetis metis(nel,nindep);
	metis.SetElementGraph(elgraph,elgraphindex);
	metis.Resequence(perm,iperm);
	fCompMesh->Permute(iperm);
*/
}

void TPZAnalysis::Assemble() {
	if(!fCompMesh || !fStructMatrix || !fSolver) 
        {
          cout << "TPZAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh << " fStructMatrix " << (void *) fStructMatrix <<
          " fSolver " << (void *) fSolver << " at file " << __FILE__ << " line " << __LINE__ << endl;
          return;
        }
        int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
        if(fSolver->Matrix() && fSolver->Matrix()->Rows()==sz)
        {
          fSolver->Matrix()->Zero();
          fStructMatrix->Assemble(*(fSolver->Matrix()),fRhs);
        }
        else
        {
          TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs);
          //mat->Print("Rigidez");
          fSolver->SetMatrix(mat);
          //aqui TPZFMatrix n� �nula
        }
        fSolver->UpdateFrom(fSolver->Matrix());
        //fRhs.Print("Rhs");
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


void TPZAnalysis::PostProcess(TPZVec<REAL> &, std::ostream &out ){

  int neq = fCompMesh->NEquations();
  TPZVec<REAL> ux((int) neq);
  TPZVec<REAL> sigx((int) neq);
  TPZManVector<REAL,3> values(3,0.);
  fCompMesh->LoadSolution(fSolution);
  //	SetExact(&Exact);
  TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
  TPZManVector<REAL,3> errors(3);
  errors.Fill(0.0);
  int nel = elvec.NElements();
  for(int i=0;i<nel;i++) {
    TPZCompEl *el = (TPZCompEl *) elvec[i];
    if(el) {
      errors.Fill(0.0);
      el->EvaluateError(fExact, errors, 0);
      int nerrors = errors.NElements();
      values.Resize(nerrors, 0.);
      for(int ier = 0; ier < nerrors; ier++)
	values[ier] += errors[ier] * errors[ier];
    }
  }

int nerrors = errors.NElements();

if (nerrors < 3) {
  PZError << endl << "TPZAnalysis::PostProcess - At least 3 norms are expected." << endl;
  out<<endl<<"############"<<endl;
  for(int ier = 0; ier < nerrors; ier++)
    out << endl << "error " << ier << "  = " << sqrt(values[ier]);
}
else{
  out << endl << "############" << endl;
  out << endl << "true_error (Norma H1) = "  << sqrt(values[0]) << endl;
  out << endl << "L2_error (Norma L2) = "    << sqrt(values[1]) << endl;
  out << endl << "estimate (Semi-norma H1) = "    << sqrt(values[2])  <<endl;
  for(int ier = 3; ier < nerrors; ier++)
    out << endl << "other norms = " << sqrt(values[ier]) << endl;
}

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

void TPZAnalysis::DefineGraphMesh(int dim, TPZVec<char *> &scalnames, TPZVec<char *> &vecnames, const char *plotfile) {

  int dim1 = dim-1;
  TPZMaterial *mat;
  int nmat,imat;
  if(!fCompMesh)
  {
    cout << "TPZAnalysis::DefineGraphMesh fCompMesh is zero\n";
    return;
  }
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

TPZMatrixSolver *TPZAnalysis::BuildPreconditioner(EPrecond preconditioner, bool overlap)
{
  if(!fSolver || !fSolver->Matrix())
  {
#ifndef BORLAND
    cout << __FUNCTION__ << " called with uninitialized stiffness matrix\n";
#else
    cout << "TPZMatrixSolver *TPZAnalysis::BuildPreconditioner" << " called with uninitialized stiffness matrix\n";
#endif
    
  }
  if(preconditioner == EJacobi)
  {
  }
  else
  {
    TPZNodesetCompute nodeset;
    TPZStack<int> elementgraph,elementgraphindex;
//    fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
    int nindep = fCompMesh->NIndependentConnects();
    int neq = fCompMesh->NEquations();
    fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
    int nel = elementgraphindex.NElements()-1;
    TPZMetis renum(nel,nindep);
    //nodeset.Print(file,elementgraphindex,elementgraph);
    renum.ConvertGraph(elementgraph,elementgraphindex,nodeset.Nodegraph(),nodeset.Nodegraphindex());
  //   cout << "nodegraphindex " << nodeset.Nodegraphindex() << endl;
  //   cout << "nodegraph " << nodeset.Nodegraph() << endl;
    nodeset.AnalyseGraph();
    //nodeset.Print(file);
    TPZStack<int> blockgraph,blockgraphindex;
    switch(preconditioner)    
    {
      case EJacobi:
        return 0;
      case EBlockJacobi:
        nodeset.BuildNodeGraph(blockgraph,blockgraphindex);
        break;
      case  EElement:
         nodeset.BuildElementGraph(blockgraph,blockgraphindex);
         break;
      case ENodeCentered:
         nodeset.BuildVertexGraph(blockgraph,blockgraphindex);
         break;
    }
    TPZStack<int> expblockgraph,expblockgraphindex;
    
    nodeset.ExpandGraph(blockgraph,blockgraphindex,fCompMesh->Block(),expblockgraph,expblockgraphindex);
    if(overlap && !(preconditioner == EBlockJacobi))
    {
      TPZSparseBlockDiagonal *sp = new TPZSparseBlockDiagonal(expblockgraph,expblockgraphindex,neq);
      TPZStepSolver *step = new TPZStepSolver(sp);
      step->SetDirect(ELU);
      step->SetReferenceMatrix(fSolver->Matrix());
      return step;
    }
    else if (overlap)
    {
      TPZBlockDiagonalStructMatrix blstr(fCompMesh);
      TPZBlockDiagonal *sp = new TPZBlockDiagonal();
      blstr.AssembleBlockDiagonal(*sp);
//      std::ofstream out("Direct assembly");
//      sp->Print("Directly assembled",out);
/*      int numbl = sp->NumberofBlocks();
      int ib,i,j;
      for(ib=numbl-1; ib>=0; ib--)
      {
        for(i=0; i<3; i++)
        {
          for(j=0; j<3; j++)
          {
            if(i!=j) (*sp)(3*ib+i,3*ib+j) = 0.;
          }
        }
      }
*/
      TPZStepSolver *step = new TPZStepSolver(sp);
      step->SetDirect(ELU);
      return step;
    }
    else
    {
      TPZVec<int> blockcolor;
      int numcolors = nodeset.ColorGraph(expblockgraph,expblockgraphindex,neq,blockcolor);
      return BuildSequenceSolver(expblockgraph,expblockgraphindex,neq,numcolors,blockcolor);
    }
  }
  return 0;
}

 /**
 * Build a sequence solver based on the block graph and its colors
 */
TPZMatrixSolver *TPZAnalysis::BuildSequenceSolver(TPZVec<int> &graph, TPZVec<int> &graphindex, int neq, int numcolors, TPZVec<int> &colors)
{
//  std::ofstream out("sequence.txt");
  TPZVec<TPZMatrix *> blmat(numcolors);
  TPZVec<TPZStepSolver *> steps(numcolors);
  int c;
  for(c=0; c<numcolors; c++)
  {
    blmat[c] = new TPZSparseBlockDiagonal(graph,graphindex, neq, c, colors);
//    blmat[c]->Print("Sparseblock matrix");
    steps[c] = new TPZStepSolver(blmat[c]);
    steps[c]->SetDirect(ELU);
    steps[c]->SetReferenceMatrix(fSolver->Matrix());
  }
  if(numcolors == 1) return steps[0];
  TPZSequenceSolver *result = new TPZSequenceSolver;
  result->ShareMatrix(*fSolver);
  for(c=numcolors-1; c>=0; c--)
  {
    result->AppendSolver(*steps[c]);
  }
  for(c=1; c<numcolors; c++)
  {
    steps[c]->SetReferenceMatrix(0);
    result->AppendSolver(*steps[c]);
  }
  for(c=0; c<numcolors; c++)
  {
    delete steps[c];
  }
  return result;
}
  
