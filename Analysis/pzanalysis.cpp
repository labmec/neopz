/**
 * \file
 * @brief Contains implementations of the TPZAnalysis methods.
 */

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

#include "TPZLagrangeMultiplier.h"
#include "pzstrmatrix.h"

#include "tpznodesetcompute.h"
#include "tpzsparseblockdiagonal.h"
#include "pzseqsolver.h"
#include "pzbdstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSloanRenumbering.h"
#include "TPZCutHillMcKee.h"

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
#define RENUMBER TPZSloanRenumbering()
#else
/**
 * @brief To renumbering will use sloan library.
 * @ingroup analysis
 */
#define RENUMBER TPZSloanRenumbering()
//#define RENUMBER TPZCutHillMcKee()
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


TPZAnalysis::TPZAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth, std::ostream &out) :
fGeoMesh(0), fCompMesh(0), fRhs(), fSolution(), fSolver(0), fStep(0), fTime(0.), fStructMatrix(0), fRenumber(new RENUMBER), fGuiInterface(NULL),  fTable()
{
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
    this->SetCompMesh(mesh, mustOptimizeBandwidth);
}

TPZAnalysis::TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth, std::ostream &out) :
fGeoMesh(0), fCompMesh(0), fRhs(), fSolution(), fSolver(0), fStep(0), fTime(0.), fStructMatrix(0), fRenumber(new RENUMBER), fGuiInterface(NULL),  fTable()
{
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
	this->SetCompMesh(mesh.operator ->(), mustOptimizeBandwidth);
}


void TPZAnalysis::SetCompMesh(TPZCompMesh * mesh, bool mustOptimizeBandwidth) {
    if(mesh)
    {
        fCompMesh = mesh;
        fGeoMesh = mesh->Reference();
        fGraphMesh[0] = 0;
        fGraphMesh[1] = 0;
        fGraphMesh[2] = 0;
        if(fSolver) fSolver->ResetMatrix();
        fCompMesh->ExpandSolution();
        long neq = fCompMesh->NEquations();
        if(neq > 20000)
        {
            std::cout << __PRETTY_FUNCTION__ << " optimizing bandwidth\n";
            std::cout.flush();
        }
        if(mustOptimizeBandwidth)
        {
            OptimizeBandwidth();
        }
        if(neq > 20000)
        {
            std::cout << __PRETTY_FUNCTION__ << " optimizing bandwidth finished\n";
            std::cout.flush();
        }
        fSolution = fCompMesh->Solution();
        fSolution.Resize(neq,1);
    }
    else
    {
        CleanUp();
    }
    fStep = 0;
    fTime = 0.;
    if(!this->fSolver){
        //seta default do stepsolver como LU
        TPZStepSolver<STATE> defaultSolver;
        defaultSolver.SetDirect(ELU);
        this->SetSolver(defaultSolver);
      
    }
    if(!this->fStructMatrix && mesh)
    {
        //seta default do StructMatrix como Full Matrix
        TPZSkylineNSymStructMatrix  defaultMatrix(mesh);
        this->SetStructuralMatrix(defaultMatrix);
    }
  

}

TPZAnalysis::~TPZAnalysis(void){
    CleanUp();
}

/// deletes all data structures
void TPZAnalysis::CleanUp()
{
    if(fSolver) delete fSolver;
    int dim;
    for(dim=0; dim<3; dim++) {
        if(fGraphMesh[dim]) delete fGraphMesh[dim];
        fGraphMesh[dim] = 0;
        fScalarNames[dim].resize(0);
        fVectorNames[dim].resize(0);
    }
    fCompMesh = 0;
    fGeoMesh = 0;
    fSolution.Redim(0,0);
    fRhs.Redim(0,0);
    fStructMatrix = 0;
    fRenumber = 0;
    fGuiInterface = 0;
    
}


void TPZAnalysis::OptimizeBandwidth() {
	//enquanto nao compilamos o BOOST no windows, vai o sloan antigo
#ifdef WIN32
	if(!fCompMesh) return;
	fCompMesh->InitializeBlock();
	TPZVec<long> perm,iperm;
	
	TPZStack<long> elgraph;
	TPZStack<long> elgraphindex;
	long nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	long nel = elgraphindex.NElements()-1;
	TPZSloan sloan(nel,nindep);
	sloan.SetElementGraph(elgraph,elgraphindex);
	sloan.Resequence(perm,iperm);
	fCompMesh->Permute(perm);
#else
	if(!fCompMesh) return;
	fCompMesh->InitializeBlock();
	
	TPZVec<long> perm,iperm;
	
	TPZStack<long> elgraph,elgraphindex;
	long nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	long nel = elgraphindex.NElements()-1;
	long el,ncel = fCompMesh->NElements();
	int maxelcon = 0;
	for(el = 0; el<ncel; el++)
	{
		TPZCompEl *cel = fCompMesh->ElementVec()[el];
		if(!cel) continue;
		std::set<long> indepconlist,depconlist;
		cel->BuildConnectList(indepconlist,depconlist);
		long locnindep = indepconlist.size();
		maxelcon = maxelcon < locnindep ? locnindep : maxelcon;
	}
	fRenumber->SetElementsNodes(nel,nindep);
	fRenumber->SetElementGraph(elgraph,elgraphindex);
	fRenumber->Resequence(perm,iperm);
	fCompMesh->Permute(perm);
    if (nel > 100000) {
        std::cout << "Applying Saddle Permute\n";
    }
    fCompMesh->SaddlePermute();
//    fCompMesh->SaddlePermute2();
	
#endif
	
}

/** @brief Determine the number of load cases from the material objects and return its value */
/**
 * this method will modify the material objects so that they have all the same number of load cases
 * the number of load cases is the maximum value of load cases of all material objects
 */
int TPZAnalysis::ComputeNumberofLoadCases()
{
    int res = 1;
    if(!fCompMesh) 
    {
        return res;
    }
    std::map<int, TPZMaterial *>::iterator it;
    // compute the maximum number of load cases for all material objects
    for( it = fCompMesh->MaterialVec().begin(); it != fCompMesh->MaterialVec().end(); it++)
    {
        TPZMaterial *mat = it->second;
        int matnumstate = mat->MinimumNumberofLoadCases();
        res = res < matnumstate ? matnumstate : res;
    }
    // set the number of load cases for all material objects
    for( it = fCompMesh->MaterialVec().begin(); it != fCompMesh->MaterialVec().end(); it++)
    {
        TPZMaterial *mat = it->second;
        mat->SetNumLoadCases(res);
    }
    return res;
}


void TPZAnalysis::AssembleResidual(){
    int numloadcases = ComputeNumberofLoadCases();

    if(Mesh()->NEquations()==0)
    {
        cout << "\n ########################################################################" <<endl;
        cout << "\n Imprimindo malha computacional de pos processamento no assemble residual " <<endl;
        cout << "\n Malha nula! " <<endl;
        DebugStop();
    }
	long sz = this->Mesh()->NEquations();
	this->Rhs().Redim(sz,numloadcases);
	fStructMatrix->Assemble(this->Rhs(),fGuiInterface);
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
    int numloadcases = ComputeNumberofLoadCases();
	long sz = fCompMesh->NEquations();
	fRhs.Redim(sz,numloadcases);
	if(fSolver->Matrix() && fSolver->Matrix()->Rows()==sz)
	{
		fSolver->Matrix()->Zero();
		fStructMatrix->Assemble(*(fSolver->Matrix().operator ->()),fRhs,fGuiInterface);
	}
	else
	{
        
		TPZMatrix<STATE> *mat = fStructMatrix->CreateAssemble(fRhs,fGuiInterface);
		fSolver->SetMatrix(mat);
		//aqui TPZFMatrix<STATE> nao eh nula
	}
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        PrintVectorByElement(sout, fRhs, 1.e-6);
//        fRhs.Print("Rhs",sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
		
	fSolver->UpdateFrom(fSolver->Matrix());
}

void TPZAnalysis::Solve() {
	long numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) 
    {
        DebugStop();
    }
	long nReducedEq = fStructMatrix->NReducedEquations();
    if (nReducedEq == numeq) 
    {
        TPZFMatrix<STATE> residual(fRhs);
        TPZFMatrix<STATE> delu(numeq,1,0.);
        //      STATE normres  = Norm(residual);
        //	cout << "TPZAnalysis::Solve residual : " << normres << " neq " << numeq << endl;
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            TPZFMatrix<STATE> res2(fRhs);
            fSolver->Matrix()->Residual(fSolution,fRhs,res2);
            std::stringstream sout;
            sout << "Residual norm " << Norm(res2) << std::endl;
    //		res2.Print("Residual",sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
    
//        {
//            std::ofstream out("Matrix.nb");
//            fSolver->Matrix()->Print("Stiffness = ",out,EMathematicaInput);
//    
//        }
        fSolver->Solve(residual, delu);
        fSolution = delu;
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            if(!fSolver->Matrix()->IsDecomposed())
            {
                TPZFMatrix<STATE> res2(fRhs);
                fSolver->Matrix()->Residual(delu,fRhs,res2);
                std::stringstream sout;
                sout << "Residual norm " << Norm(res2) << std::endl;
                //            res2.Print("Residual",sout);
                LOGPZ_DEBUG(logger,sout.str())
            }
        }
#endif
    
    }
    else 
    {
        TPZFMatrix<STATE> residual(nReducedEq,1,0.);
    	TPZFMatrix<STATE> delu(nReducedEq,1,0.);
        fStructMatrix->EquationFilter().Gather(fRhs,residual);
	    fSolver->Solve(residual, delu);
        fSolution.Redim(numeq,1);
        fStructMatrix->EquationFilter().Scatter(delu,fSolution);
    }
#ifdef LOG4CXX
    std::stringstream sout;
    TPZStepSolver<STATE> *step = dynamic_cast<TPZStepSolver<STATE> *> (fSolver);
    if(!step) DebugStop();
    long nsing = step->Singular().size();
	if(nsing && logger->isWarnEnabled()) {
		sout << "Number of singular equations " << nsing;
		std::list<long>::iterator it = step->Singular().begin();
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
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Solution norm " << Norm(fSolution) << std::endl;
		fSolution.Print("delu",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif	
	fCompMesh->LoadSolution(fSolution);
    fCompMesh->TransferMultiphysicsSolution();

}

void TPZAnalysis::LoadSolution() {	
	if(fCompMesh) {
		fCompMesh->LoadSolution(fSolution);
	}
}

void TPZAnalysis::Print( const std::string &name, std::ostream &out) {
	out<<endl<<name<<endl;
    if (!fCompMesh) {
        out << "No computational mesh\n";
    }
    else
    {
        long i,nelements = fCompMesh->ConnectVec().NElements();
        for(i=0;i<nelements;i++) {
            TPZConnect &gnod = fCompMesh->ConnectVec()[i];
            if(gnod.SequenceNumber()!=-1) {
                out << "Connect index " << i << endl;
                gnod.Print(*fCompMesh,out);
            }
        }
    }
	fSolution.Print("fSolution",out);
    if (fCompMesh) {
        fCompMesh->ConnectSolution(out);
    }
}


void TPZAnalysis::PostProcess(TPZVec<REAL> &ervec, std::ostream &out) {
	long i;
	//TPZVec<REAL> ux((int) neq);
	//TPZVec<REAL> sigx((int) neq);
	TPZManVector<REAL,10> values(10,0.);
	TPZManVector<REAL,10> values2(10,0.);
	fCompMesh->LoadSolution(fSolution);
	TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
	TPZManVector<REAL,10> errors(10);
	errors.Fill(0.0);
	long nel = elvec.NElements();
	int matId0 = 0;
	for(i=0;i<nel;i++) {
        TPZCompEl *cel = elvec[i];
        if(!cel) continue;
		matId0=cel->Material()->Id();
		if(matId0 > 0)
			break;
	}
	bool lastEl=false;
	for(i=0;i<nel;i++) {
		TPZCompEl *el = (TPZCompEl *) elvec[i];
		if(el) {
			errors.Fill(0.0);
            TPZGeoEl *gel = el->Reference();
            if(gel->Dimension() != fCompMesh->Dimension()) continue;
			el->EvaluateError(fExact, errors, 0);
			if(matId0==el->Material()->Id()){
				for(int ier = 0; ier < errors.NElements(); ier++) 	values[ier] += errors[ier] * errors[ier];
				lastEl=false;
			}
			else {
				for(int ier = 0; ier < errors.NElements(); ier++)	values2[ier] += errors[ier] * errors[ier];
				lastEl=true;
			}
		}
	}
	int nerrors = errors.NElements();
	ervec.Resize(2*nerrors);
	ervec.Fill(-10.0);
	
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
		
		for(i=0;i<nerrors;i++) {
			ervec[i] = sqrt(values[i]);
		}
		for(i=nerrors;i<2*nerrors;i++){
			ervec[i] = sqrt(values2[i-nerrors]);
		}
	}
	return;
}

void TPZAnalysis::PostProcessError(TPZVec<REAL> &ervec, std::ostream &out ){
    
    long neq = fCompMesh->NEquations();
    TPZVec<REAL> ux(neq);
    TPZVec<REAL> sigx(neq);
    TPZManVector<REAL,10> values(10,0.);
    fCompMesh->LoadSolution(fSolution);
    //	SetExact(&Exact);
    TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
    TPZManVector<REAL,10> errors(10);
    errors.Fill(0.0);
    long i, nel = elvec.NElements();
    for(i=0;i<nel;i++) {
        TPZCompEl *el = (TPZCompEl *) elvec[i];
        if(el) {
            TPZMaterial *mat = el->Material();
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if(!bc)
            {
                errors.Fill(0.0);
                el->EvaluateError(fExact, errors, 0);
                int nerrors = errors.NElements();
                values.Resize(nerrors, 0.);
                for(int ier = 0; ier < nerrors; ier++)
                {
                    values[ier] += errors[ier] * errors[ier];
                }
            }
        }
    }
    
    int nerrors = errors.NElements();
	ervec.Resize(nerrors);
	ervec.Fill(-10.0);

    if (nerrors < 3) {
        PZError << endl << "TPZAnalysis::PostProcess - At least 3 norms are expected." << endl;
        out<<endl<<"############"<<endl;
        for(int ier = 0; ier < nerrors; ier++)
            out << endl << "error " << ier << "  = " << sqrt(values[ier]);
    }
    else{
        out << "############" << endl;
        out <<"Norma H1 = "  << sqrt(values[0]) << endl;
        out <<"Norma L2 = "    << sqrt(values[1]) << endl;
        out << "Semi-norma H1 = "    << sqrt(values[2])  <<endl;
        for(int ier = 3; ier < nerrors; ier++)
            out << "other norms = " << sqrt(values[ier]) << endl;
    }
	// Returns the calculated errors.
	for(i=0;i<nerrors;i++)
		ervec[i] = sqrt(values[i]);
    return;
}

void TPZAnalysis::PostProcessTable( TPZFMatrix<REAL> &,std::ostream & )//pos,out
{
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	long nel = elvec.NElements();
	for(long i=0;i<nel;i++) {
		//TPZCompEl *el = (TPZCompEl *) elvec[i];
		//if(el)	el->PrintTable(fExact,pos,out);
	}
	return;
}

void TPZAnalysis::ShowShape(const std::string &plotfile, TPZVec<long> &equationindices)
{
	
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("Solution");
    DefineGraphMesh(fCompMesh->Dimension(), scalnames, vecnames, plotfile);
    int porder = fCompMesh->GetDefaultOrder();
    
    int neq = equationindices.size();
    TPZFMatrix<STATE> solkeep(fSolution);
    fSolution.Zero();
    for (int ieq = 0; ieq < neq; ieq++) {
        fSolution(equationindices[ieq],0) = 1.;
        LoadSolution();
        Mesh()->TransferMultiphysicsSolution();
        PostProcess(porder+1);
        fSolution.Zero();
    }
    fSolution = solkeep;
    LoadSolution();
}

void TPZAnalysis::LoadShape(double ,double , long ,TPZConnect* start){
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
    long neq = fCompMesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
	Assemble();

    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }

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
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		TPZBndCond *bc = dynamic_cast<TPZBndCond *> (matit->second);
        TPZMaterial *mat = matit->second;
		if(mat && !bc && mat->Dimension() == dim) break;
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
	int posplot = plotfile.rfind(".plt");
	int posdx = plotfile.rfind(".dx");
	int pospos = plotfile.rfind(".pos");
	int posvtk = plotfile.rfind(".vtk");
	long filelength = plotfile.size();
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
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		TPZBndCond *bc = dynamic_cast<TPZBndCond *>(matit->second);
        TPZLagrangeMultiplier *lag = dynamic_cast<TPZLagrangeMultiplier *>(matit->second);
		if(matit->second && !bc && !lag && matit->second->Dimension() == dimension) break;
	}
	if(matit == fCompMesh->MaterialVec().end()) return;
	fGraphMesh[dim1]->SetCompMesh(fCompMesh,matit->second);
	
	fGraphMesh[dim1]->SetResolution(resolution);
	fGraphMesh[dim1]->DrawMesh(1);
	fGraphMesh[dim1]->DrawSolution(fStep,fTime);
	fStep++;
}

void TPZAnalysis::AnimateRun(long num_iter, int steps, TPZVec<std::string> &scalnames,
							 TPZVec<std::string> &vecnames, const std::string &plotfile) {
	Assemble();
	
	long numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return;
	
	TPZFMatrix<STATE> residual(fRhs);
	int dim = HighestDimension();
	TPZMaterial * mat = 0;
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		TPZBndCond *bc = dynamic_cast<TPZBndCond *>(matit->second);
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
	gg.SetFileName(plotfile);
	gg.SetResolution(0);
	gg.DrawMesh(num_iter);
	
	long i;
	for(i=1; i<=num_iter;i+=steps){
		
		
		TPZStepSolver<STATE> sol;
		sol.ShareMatrix(Solver());
		sol.SetJacobi(i,0.,0);
		SetSolver(sol);
		fSolver->Solve(fRhs, fSolution);
		
		fCompMesh->LoadSolution(fSolution);
		gg.DrawSolution(i-1,0);
	}
}

int TPZAnalysis::HighestDimension(){
	int dim = 0;
	std::map<int, TPZMaterial * >::iterator matit;
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

void TPZAnalysis::DefineElementTable(int dimension, TPZVec<long> &GeoElIds, TPZVec<REAL> &points) {
	fTable.fDimension = dimension;
	fTable.fGeoElId = GeoElIds;
	fTable.fLocations = points;
	long numel = GeoElIds.NElements();
	fTable.fCompElPtr.Resize(numel);
	long iel;
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
		long numel = fTable.fCompElPtr.NElements();
		for(long iel=0; iel<numel; iel++) {
			cel = (TPZCompEl *) fTable.fCompElPtr[iel];
			if(cel) cel->PrintTitle((char *)fTable.fVariableNames[iv],*(fTable.fOutfile));
		}
	}
	*(fTable.fOutfile) << endl;
	int dim;
	TPZVec<REAL> point(fTable.fDimension);
	for(dim=1; dim<fTable.fDimension+1; dim++) {
		for(int iv=0; iv<numvar; iv++) {
			long numel = fTable.fCompElPtr.NElements();
			for(long iel=0; iel<numel; iel++) {
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
		long numel = fTable.fCompElPtr.NElements();
		for(long iel=0; iel<numel; iel++) {
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
void TPZAnalysis::SetSolver(TPZMatrixSolver<STATE> &solver){
	if(fSolver) delete fSolver;
    fSolver = (TPZMatrixSolver<STATE> *) solver.Clone();
}

TPZMatrixSolver<STATE> *TPZAnalysis::BuildPreconditioner(EPrecond preconditioner, bool overlap)
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
		TPZStack<long> elementgraph,elementgraphindex;
		long nindep = fCompMesh->NIndependentConnects();
		long neq = fCompMesh->NEquations();
		fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
		long nel = elementgraphindex.NElements()-1;
		TPZMetis renum(nel,nindep);
		renum.ConvertGraph(elementgraph,elementgraphindex,nodeset.Nodegraph(),nodeset.Nodegraphindex());
		nodeset.AnalyseGraph();

		TPZStack<long> blockgraph,blockgraphindex;
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
		TPZStack<long> expblockgraph,expblockgraphindex;
		
		nodeset.ExpandGraph(blockgraph,blockgraphindex,fCompMesh->Block(),expblockgraph,expblockgraphindex);
#ifdef LOG4CXX
#ifdef PZDEBUG2
        if (logger->isDebugEnabled())
        {
            std::map<long,long> blocksizes;
            long i;
            long totalsize = 0;
            for(i=0; i< expblockgraphindex.NElements()-1;i++)
            {
                long bls = expblockgraphindex[i+1]-expblockgraphindex[i];
                blocksizes[bls]++;
                totalsize += bls*bls;
            }
            std::map<long,long>::iterator it;
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " total size of allocation " << totalsize << std::endl;
            for(it=blocksizes.begin(); it != blocksizes.end(); it++)
            {
                sout << "block size " << (*it).first << " number of blocks " << (*it).second << std::endl;
            }
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
#endif
#endif
		if(overlap && !(preconditioner == EBlockJacobi))
		{
			TPZSparseBlockDiagonal<STATE> *sp = new TPZSparseBlockDiagonal<STATE>(expblockgraph,expblockgraphindex,neq);
			TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(sp);
			step->SetDirect(ELU);
			step->SetReferenceMatrix(fSolver->Matrix());
			return step;
		}
		else if (overlap)
		{
			TPZBlockDiagonalStructMatrix blstr(fCompMesh);
			TPZBlockDiagonal<STATE> *sp = new TPZBlockDiagonal<STATE>();
			blstr.AssembleBlockDiagonal(*sp);
			TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(sp);
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

/** @brief Build a sequence solver based on the block graph and its colors */
TPZMatrixSolver<STATE> *TPZAnalysis::BuildSequenceSolver(TPZVec<long> &graph, TPZVec<long> &graphindex, long neq, int numcolors, TPZVec<int> &colors)
{
	TPZVec<TPZMatrix<STATE> *> blmat(numcolors);
	TPZVec<TPZStepSolver<STATE> *> steps(numcolors);
	int c;
	for(c=0; c<numcolors; c++)
	{
		blmat[c] = new TPZSparseBlockDiagonal<STATE>(graph,graphindex, neq, c, colors);
		steps[c] = new TPZStepSolver<STATE>(blmat[c]);
		steps[c]->SetDirect(ELU);
		steps[c]->SetReferenceMatrix(fSolver->Matrix());
	}
	if(numcolors == 1) return steps[0];
	TPZSequenceSolver<STATE> *result = new TPZSequenceSolver<STATE>;
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

/// Integrate the postprocessed variable name over the elements included in the set matids
TPZVec<STATE> TPZAnalysis::Integrate(const std::string &varname, const std::set<int> &matids)
{
    // the postprocessed index of the varname for each material id
    std::map<int,int> variableids;
    int nvars = 0;
    
    std::map<int,TPZMaterial *>::iterator itmap;
    for (itmap = fCompMesh->MaterialVec().begin(); itmap != fCompMesh->MaterialVec().end(); itmap++) {
        if (matids.find(itmap->first) != matids.end()) {
            variableids[itmap->first] = itmap->second->VariableIndex(varname);
            nvars = itmap->second->NSolutionVariables(variableids[itmap->first]);
        }
    }
    TPZManVector<STATE,3> result(nvars,0.);
    long nelem = fCompMesh->NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCompMesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int matid = gel->MaterialId();
        if (matids.find(matid) == matids.end()) {
            continue;
        }
        TPZManVector<STATE,3> locres(nvars,0.);
        locres = cel->IntegrateSolution(variableids[matid]);
        for (int iv=0; iv<nvars; iv++)
        {
            result[iv] += locres[iv];
        }
    }
    return result;
}

/// extract the values corresponding to the connect from the vector
static void ConnectSolution(long cindex, TPZCompMesh *cmesh, TPZFMatrix<STATE> &glob, TPZVec<STATE> &sol)
{
    long seqnum = cmesh->ConnectVec()[cindex].SequenceNumber();
    int blsize = cmesh->Block().Size(seqnum);
    int position = cmesh->Block().Position(seqnum);
    sol.resize(blsize);
    for (long i=position; i< position+blsize; i++) {
        sol[i-position] = glob(i,0);
    }
}

static STATE ConnectNorm(long cindex, TPZCompMesh *cmesh, TPZFMatrix<STATE> &glob)
{
    TPZManVector<STATE,20> cvec;
    ConnectSolution(cindex, cmesh, glob, cvec);
    STATE norm = 0.;
    for (long i=0; i<cvec.size(); i++) {
        norm += cvec[i]*cvec[i];
    }
    norm = sqrt(norm);
    return norm;
}

/// Print the residual vector for those elements with entry above a given tolerance
void TPZAnalysis::PrintVectorByElement(std::ostream &out, TPZFMatrix<STATE> &vec, REAL tol)
{
    long nel = fCompMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fCompMesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        int ic;
        for (ic=0; ic<nc; ic++) {
            long cindex = cel->ConnectIndex(ic);
            TPZConnect &c = cel->Connect(ic);
            if(c.HasDependency() || c.IsCondensed())
            {
                continue;
            }
#ifndef STATE_COMPLEX
            if (ConnectNorm(cindex, fCompMesh, vec) >= tol) {
                break;
            }
#endif
        }
        // if all connects have norm below tol, do not print the element
        if (ic == nc) {
            continue;
        }
        bool hasgeometry = false;
        TPZManVector<REAL> xcenter(3,0.);
        TPZGeoEl *gel = cel->Reference();
        // if a geometric element exists, extract its center coordinate
        if (gel) {
            hasgeometry = true;
            TPZManVector<REAL,3> xicenter(gel->Dimension(),0.);
            gel->CenterPoint(gel->NSides()-1, xicenter);
            gel->X(xicenter,xcenter);
        }
        out << "CompEl " << el;
        if (hasgeometry) {
            out << " Gel " << gel->Index() << " matid " << gel->MaterialId() << " Center " << xcenter << std::endl;
        }
        for (ic = 0; ic<nc; ic++) {
            TPZManVector<STATE> connectsol;
            long cindex = cel->ConnectIndex(ic);
            TPZConnect &c = fCompMesh->ConnectVec()[cindex];
            if(c.HasDependency() || c.IsCondensed()) {
                out << "connect " << ic << " is restrained\n";
                continue;
            }
            ConnectSolution(cindex, fCompMesh, vec, connectsol);
            for (int i=0; i<connectsol.size(); i++) {
                if (fabs(connectsol[i]) < tol) {
                    connectsol[i] = 0.;
                }
            }
            out << ic << " index " << cindex << " values " << connectsol << std::endl;
        }
    }

}

