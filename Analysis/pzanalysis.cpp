/**
 * \file
 * @brief Contains implementations of the TPZAnalysis methods.
 */
//$Id: pzanalysis.cpp,v 1.59 2011-05-11 01:12:37 phil Exp $
//$Id: pzanalysis.cpp,v 1.59 2011-05-11 01:12:37 phil Exp $

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
#include "pzvtkmesh.h"


#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzmetis.h"
#include "pzsloan.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzstrmatrix.h"

#include "tpznodesetcompute.h"
#include "tpzsparseblockdiagonal.h"
#include "pzseqsolver.h"
#include "pzbdstrmatrix.h"


#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.analysis"));
#endif

#ifdef USING_BOOST
#include "TPZBoostGraph.h"
/**
 * @brief To renumbering will use boost library.
 * @ingroup analysis
 */
#define RENUMBER TPZBoostGraph(TPZBoostGraph::KMCExpensive)
#else
/**
 * @brief To renumbering will use sloan library.
 * @ingroup analysis
 */
#define RENUMBER TPZSloan()
#endif

#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void TPZAnalysis::SetStructuralMatrix(TPZStructMatrix &strmatrix){
	fStructMatrix = TPZAutoPointer<TPZStructMatrix>(strmatrix.Clone());
}

void TPZAnalysis::SetStructuralMatrix(TPZAutoPointer<TPZStructMatrix> strmatrix){
	fStructMatrix = TPZAutoPointer<TPZStructMatrix>(strmatrix->Clone());
}
TPZAnalysis::TPZAnalysis() : fGeoMesh(0), fCompMesh(0), fRhs(), fSolution(), fSolver(0), fStep(0), fTime(0.), fStructMatrix(0), fRenumber(new RENUMBER)
, fGuiInterface(NULL), fTable() {
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
}


TPZAnalysis::TPZAnalysis(TPZCompMesh *mesh,std::ostream &out) :
fGeoMesh(0), fCompMesh(0), fRhs(), fSolution(), fSolver(0), fStep(0), fTime(0.), fStructMatrix(0), fRenumber(new RENUMBER), fGuiInterface(NULL),  fTable()
{
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
	this->SetCompMesh(mesh);
}

TPZAnalysis::TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh,std::ostream &out) :
fGeoMesh(0), fCompMesh(0), fRhs(), fSolution(), fSolver(0), fStep(0), fTime(0.), fStructMatrix(0), fRenumber(new RENUMBER), fGuiInterface(NULL),  fTable()
{
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
	this->SetCompMesh(mesh.operator ->());
}


void TPZAnalysis::SetCompMesh(TPZCompMesh * mesh) {
	fCompMesh = mesh;
	fGeoMesh = mesh->Reference();
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
	if(fSolver) fSolver->ResetMatrix();
	SetBlockNumber();
	fStep = 0;
	fTime = 0.;
	fCompMesh->ExpandSolution();
	fSolution = fCompMesh->Solution();
	int neq = fCompMesh->NEquations();
	fSolution.Resize(neq,1);
}

TPZAnalysis::~TPZAnalysis(void){
	if(fSolver) delete fSolver;
	int dim;
	for(dim=0; dim<3; dim++) {
		if(fGraphMesh[dim]) delete fGraphMesh[dim];
		fGraphMesh[dim] = 0;
	}
}

void TPZAnalysis::SetBlockNumber(){
	//enquanto nao compilamos o BOOST no windows, vai o sloan antigo
#ifdef WIN32
	if(!fCompMesh) return;
	fCompMesh->InitializeBlock();
	TPZVec<int> perm,iperm;
	
	TPZStack<int> elgraph,elgraphindex;
	int nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int nel = elgraphindex.NElements()-1;
	TPZSloan sloan(nel,nindep);
	sloan.SetElementGraph(elgraph,elgraphindex);
	sloan.Resequence(perm,iperm);
	fCompMesh->Permute(perm);
#else
	if(!fCompMesh) return;
	fCompMesh->InitializeBlock();
	
	TPZVec<int> perm,iperm;
	
	TPZStack<int> elgraph,elgraphindex;
	int nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int nel = elgraphindex.NElements()-1;
	int el,ncel = fCompMesh->NElements();
	int maxelcon = 0;
	for(el = 0; el<ncel; el++)
	{
		TPZCompEl *cel = fCompMesh->ElementVec()[el];
		if(!cel) continue;
		std::set<int> indepconlist,depconlist;
		cel->BuildConnectList(indepconlist,depconlist);
		int locnindep = indepconlist.size();
		maxelcon = maxelcon < locnindep ? locnindep : maxelcon;
	}
	fRenumber->SetElementsNodes(nel,nindep);
	//	TPZSloan sloan(nel,nindep,maxelcon);
	fRenumber->SetElementGraph(elgraph,elgraphindex);
	fRenumber->Resequence(perm,iperm);
	fCompMesh->Permute(perm);
	/*
	 fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	 
	 TPZMetis metis(nel,nindep);
	 metis.SetElementGraph(elgraph,elgraphindex);
	 metis.Resequence(perm,iperm);
	 fCompMesh->Permute(iperm);
	 */
	
#endif
	
}

void TPZAnalysis::AssembleResidual(){
	int sz = this->Mesh()->NEquations();
	this->Rhs().Redim(sz,1);
	fStructMatrix->Assemble(this->Rhs(),fGuiInterface);
	//TPZStructMatrix::Assemble(this->Rhs(), *this->Mesh());
}//void

void TPZAnalysis::Assemble()
{
	if(!fCompMesh || !fStructMatrix || !fSolver)
	{
		std::stringstream sout;
		sout << "TPZAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh
		<< " fStructMatrix " << (void *) fStructMatrix.operator->()
		<< " fSolver " << (void *) fSolver;
#ifndef WINDOWS
		sout << " at file " << __FILE__ << " line " << __LINE__ ;
#else
		sout << " TPZAnalysis::Assemble() " ;
#endif
#ifdef LOG4CXX
		LOGPZ_ERROR(logger,sout.str().c_str());
#else
		std::cout << sout.str().c_str() << std::endl;
#endif
		return;
	}
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	if(fSolver->Matrix() && fSolver->Matrix()->Rows()==sz)
	{
		fSolver->Matrix()->Zero();
		fStructMatrix->Assemble(*(fSolver->Matrix().operator ->()),fRhs,fGuiInterface);
	}
	else
	{
		TPZMatrix<REAL> *mat = fStructMatrix->CreateAssemble(fRhs,fGuiInterface);
		fSolver->SetMatrix(mat);
		//aqui TPZFMatrix<REAL> n� �nula
	}
	
	//   ofstream fileout("rigidez.txt");
	//   fSolver->Matrix()->Print("Rigidez", fileout, EMathematicaInput);  
	
	
	fSolver->UpdateFrom(fSolver->Matrix());
	
	//fRhs.Print("Rhs");
	//cout.flush();
}

void TPZAnalysis::Solve() {
	int numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return;
	
	TPZFMatrix<REAL> residual(fRhs);
	TPZFMatrix<REAL> delu(numeq,1,0.);
	/*	if(fSolution.Rows() != numeq) {
	 fSolution.Redim(numeq,1);
	 } else {
	 fSolver->Matrix()->Residual(fSolution,fRhs,residual);
	 }*/
	//      REAL normres  = Norm(residual);
	//	cout << "TPZAnalysis::Solve residual : " << normres << " neq " << numeq << endl;
#ifdef LOG4CXX_KEEP
	{
		std::stringstream sout;
		sout << "Residual norm " << Norm(residual) << std::endl;
		residual.Print("Residual",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fSolver->Solve(residual, delu);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		TPZStepSolver<REAL> *step = dynamic_cast<TPZStepSolver<REAL> *> (fSolver);
		if(!step) DebugStop();
		int nsing = step->Singular().size();
		sout << "Number of singular equations " << nsing;
		std::list<int>::iterator it = step->Singular().begin();
		if(nsing) sout << "\nSingular modes ";
		while(it != step->Singular().end())
		{
			sout << *it << " ";
			it++;
		}
		if(nsing) sout << std::endl;
		LOGPZ_WARN(logger,sout.str())
	}
#endif
#ifdef LOG4CXX_KEEP
	{
		std::stringstream sout;
		sout << "Solution norm " << Norm(delu) << std::endl;
		delu.Print("delu",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fSolution = delu;
	
	fCompMesh->LoadSolution(fSolution);
}

void TPZAnalysis::LoadSolution(){
	
	if(fCompMesh) {
		fCompMesh->LoadSolution(fSolution);
	}
}

void TPZAnalysis::Print( const std::string &name , std::ostream &out )
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
	
	int i, neq = fCompMesh->NEquations();
	TPZVec<REAL> ux((int) neq);
	TPZVec<REAL> sigx((int) neq);
	TPZManVector<REAL,10> values(10,0.);
	TPZManVector<REAL,10> values2(10,0.);
	fCompMesh->LoadSolution(fSolution);
	//	SetExact(&Exact);
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	TPZManVector<REAL,10> errors(10);
	errors.Fill(0.0);
	int nel = elvec.NElements();
	int matId0;
	for(i=0;i<nel;i++) {
		matId0=elvec[i]->Material()->Id();
		if(matId0 > 0)
			break;
	}
	bool lastEl=false;
	for(i=0;i<nel;i++) {
		TPZCompEl *el = (TPZCompEl *) elvec[i];
		if(el) {
			errors.Fill(0.0);
			el->EvaluateError(fExact, errors, 0);
			if(matId0==el->Material()->Id()){
				for(int ier = 0; ier < errors.NElements(); ier++) 	values[ier] += errors[ier] * errors[ier];
				lastEl=false;
			}
			else{
				for(int ier = 0; ier < errors.NElements(); ier++)	values2[ier] += errors[ier] * errors[ier];
				lastEl=true;
			}
		}
		
	}
	int nerrors = errors.NElements();
	
	if (nerrors==4) {
		if(lastEl){
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values2[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values2[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values2[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values2[3])  <<endl;

			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values[2])  <<endl;
			
		}
		else{
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values[3])  <<endl;
			
			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values2[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values2[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values2[2])  <<endl;
		}
	}
	else {
		if(lastEl){
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values[3])  <<endl;
			
			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values2[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values2[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values2[2])  <<endl;
			
		}
		else{
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values2[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values2[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values2[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values2[3])  <<endl;
			
			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values[2])  <<endl;
		}
	}
	return;
}

void TPZAnalysis::PostProcessTable( TPZFMatrix<REAL> &,std::ostream & )//pos,out
{
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	int nel = elvec.NElements();
	for(int i=0;i<nel;i++) {
		//TPZCompEl *el = (TPZCompEl *) elvec[i];
		//if(el)	el->PrintTable(fExact,pos,out);
	}
	return;
}

void TPZAnalysis::ShowShape( TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames, char *plotfile, std::ostream &) {//1o :TPZConnect  *nod,
	
    TPZAutoPointer<TPZMaterial> mat = fCompMesh->MaterialVec().rbegin()->second;
	TPZV3DGraphMesh gg(fCompMesh,2,mat);
	
	gg.SetFileName(plotfile);
	gg.SetResolution(0);
	gg.DrawMesh(1);
	gg.SetNames(scalnames,vecnames);
	gg.DrawSolution(0,0.);
	
}

void TPZAnalysis::LoadShape(double ,double , int ,TPZConnect* start){
	//void TPZAnalysis::LoadShape(double dx,double dy, int numelem,TPZConnect* start){
	Assemble();
	fRhs.Zero();
	fSolution.Zero();
	
}
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void TPZAnalysis::Run(std::ostream &out)
{
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
	Assemble();
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    Solve();
	
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
	
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
}

void TPZAnalysis::DefineGraphMesh(int dim, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile) {
	
	int dim1 = dim-1;
	if(!fCompMesh)
	{
		cout << "TPZAnalysis::DefineGraphMesh fCompMesh is zero\n";
		return;
	}
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		TPZBndCond *bc = dynamic_cast<TPZBndCond *> (matit->second.operator->());
		if(matit->second && !bc && matit->second->Dimension() == dim) break;
	}
	if(matit == fCompMesh->MaterialVec().end())
	{
		std::cout << __PRETTY_FUNCTION__ << " The computational mesh has no associated material!!!!\n";
		DebugStop();
		return;
	}
	if(fGraphMesh[dim1]) delete fGraphMesh[dim1];
	fScalarNames[dim1] = scalnames;
	fVectorNames[dim1] = vecnames;
	//misael
	//  if(! ( strcmp( strrchr(plotfile,'.')+1,"plt") ) )	{
	int posplot = plotfile.rfind(".plt");
	int posdx = plotfile.rfind(".dx");
	int pospos = plotfile.rfind(".pos");
	int posvtk = plotfile.rfind(".vtk");
	int filelength = plotfile.size();
	if(filelength-posplot == 3)	{
		fGraphMesh[dim1] = new TPZV3DGraphMesh(fCompMesh,dim,matit->second) ;
	}else if(filelength-posdx == 3) {
		fGraphMesh[dim1] = new TPZDXGraphMesh(fCompMesh,dim,matit->second,scalnames,vecnames) ;
	}else if(filelength-pospos == 3) {
		fGraphMesh[dim1] = new TPZMVGraphMesh(fCompMesh,dim,matit->second);
	}
	else if(filelength-posvtk == 4) {
		fGraphMesh[dim1] = new TPZVTKGraphMesh(fCompMesh,dim,matit->second,scalnames,vecnames);
	} else {
		cout << "grafgrid was not created\n";
		fGraphMesh[dim1] = 0;
	}
	if(fGraphMesh[dim1]) {
		fGraphMesh[dim1]->SetFileName(plotfile);
	}
}

void TPZAnalysis::CloseGraphMesh(){
	for(int i = 0; i < 3; i++){
		if ( this->fGraphMesh[i] ){
			delete this->fGraphMesh[i];
			this->fGraphMesh[i] = NULL;
		}
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
	//  TPZMaterial *mat;
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		TPZBndCond *bc = dynamic_cast<TPZBndCond *>(matit->second.operator->());
		if(matit->second && !bc && matit->second->Dimension() == dimension) break;
	}
	if(matit == fCompMesh->MaterialVec().end()) return;
	fGraphMesh[dim1]->SetCompMesh(fCompMesh,matit->second);
	
	fGraphMesh[dim1]->SetResolution(resolution);
	fGraphMesh[dim1]->DrawMesh(1);
	fGraphMesh[dim1]->DrawSolution(fStep,fTime);
	//   delete fGraphMesh[dim1];
	//   fGraphMesh[dim1] = 0;
	fStep++;
}

void TPZAnalysis::AnimateRun(int num_iter, int steps,
							 TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames, const std::string &plotfile)
{
	Assemble();
	
	int numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return;
	
	TPZFMatrix<REAL> residual(fRhs);
	int dim = HighestDimension();
	TPZAutoPointer<TPZMaterial> mat = 0;
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		TPZBndCond *bc = dynamic_cast<TPZBndCond *>(matit->second.operator->());
		if(bc) continue;
		if( matit->second->Dimension() == dim)
		{
			mat = matit->second;
			break;
		}
	}
	if(!mat)
	{
		std::cout << __PRETTY_FUNCTION__ << " no material found " << std::endl;
		LOGPZ_ERROR(logger, " no material found");
		return;
	}
	TPZDXGraphMesh gg(fCompMesh,dim,mat,scalnames,vecnames) ;
	// 		TPZV3DGrafGrid gg(fCompMesh) ;
	gg.SetFileName(plotfile);
	gg.SetResolution(0);
	gg.DrawMesh(num_iter);
	
	int i;
	for(i=1; i<=num_iter;i+=steps){
		
		
		TPZStepSolver<REAL> sol;
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
	int dim = 0;
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		dim = matit->second->Dimension() > dim ? matit->second->Dimension() : dim;
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
void TPZAnalysis::SetSolver(TPZMatrixSolver<REAL> &solver){
	if(fSolver) delete fSolver;
    fSolver = (TPZMatrixSolver<REAL> *) solver.Clone();
}

TPZMatrixSolver<REAL> *TPZAnalysis::BuildPreconditioner(EPrecond preconditioner, bool overlap)
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
#ifdef LOG4CXX
#ifdef DEBUG2
		std::map<int,int> blocksizes;
		int i;
		int totalsize;
		for(i=0; i< expblockgraphindex.NElements()-1;i++)
		{
			int bls = expblockgraphindex[i+1]-expblockgraphindex[i];
			blocksizes[bls]++;
			totalsize += bls*bls;
		}
		std::map<int,int>::iterator it;
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " total size of allocation " << totalsize << std::endl;
		for(it=blocksizes.begin(); it != blocksizes.end(); it++)
		{
			sout << "block size " << (*it).first << " number of blocks " << (*it).second << std::endl;
		}
		LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
#endif
		if(overlap && !(preconditioner == EBlockJacobi))
		{
			TPZSparseBlockDiagonal<REAL> *sp = new TPZSparseBlockDiagonal<REAL>(expblockgraph,expblockgraphindex,neq);
			TPZStepSolver<REAL> *step = new TPZStepSolver<REAL>(sp);
			step->SetDirect(ELU);
			step->SetReferenceMatrix(fSolver->Matrix());
			return step;
		}
		else if (overlap)
		{
			TPZBlockDiagonalStructMatrix blstr(fCompMesh);
			TPZBlockDiagonal<REAL> *sp = new TPZBlockDiagonal<REAL>();
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
			TPZStepSolver<REAL> *step = new TPZStepSolver<REAL>(sp);
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
TPZMatrixSolver<REAL> *TPZAnalysis::BuildSequenceSolver(TPZVec<int> &graph, TPZVec<int> &graphindex, int neq, int numcolors, TPZVec<int> &colors)
{
	//  std::ofstream out("sequence.txt");
	TPZVec<TPZMatrix<REAL> *> blmat(numcolors);
	TPZVec<TPZStepSolver<REAL> *> steps(numcolors);
	int c;
	for(c=0; c<numcolors; c++)
	{
		blmat[c] = new TPZSparseBlockDiagonal<REAL>(graph,graphindex, neq, c, colors);
		//    blmat[c]->Print("Sparseblock matrix");
		steps[c] = new TPZStepSolver<REAL>(blmat[c]);
		steps[c]->SetDirect(ELU);
		steps[c]->SetReferenceMatrix(fSolver->Matrix());
	}
	if(numcolors == 1) return steps[0];
	TPZSequenceSolver<REAL> *result = new TPZSequenceSolver<REAL>;
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

