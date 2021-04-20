/**
 * \file
 * @brief Contains implementations of the TPZAnalysis methods.
 */

#include "pzanalysis.h"
#include <math.h>                          // for sqrt, fabs
#include <thread>
#include <vector>
#include <stdio.h>                         // for NULL
#include <string.h>                        // for strcpy, strlen
#include <chrono>
#ifdef MACOSX
#include <__functional_base>               // for less
#include <__tree>                          // for __tree_const_iterator, ope...
#endif
#include <list>                            // for list, __list_iterator, lis...
#include <map>                             // for __map_iterator, map, map<>...
#include <string>                          // for allocator, basic_string
#include <utility>                         // for pair
#include "TPZLagrangeMultiplier.h"         // for TPZLagrangeMultiplier
#include "TPZSkylineNSymStructMatrix.h"    // for TPZSkylineNSymStructMatrix
#include "TPZSloanRenumbering.h"           // for TPZSloanRenumbering
#include "pzadmchunk.h"                    // for TPZAdmChunkVector
#include "pzbdstrmatrix.h"                 // for TPZBlockDiagonalStructMatrix
#include "pzblock.h"                       // for TPZBlock
#include "pzblockdiag.h"                   // for TPZBlockDiagonal
#include "pzbndcond.h"                     // for TPZBndCond
#include "TPZChunkVector.h"                       // for TPZChunkVector
#include "pzcmesh.h"                       // for TPZCompMesh
#include "pzcompel.h"                      // for TPZCompEl
#include "pzconnect.h"                     // for TPZConnect
#include "pzdxmesh.h"                      // for TPZDXGraphMesh
#include "TPZEquationFilter.h"              // for TPZEquationFilter
#include "pzgeoel.h"                       // for TPZGeoEl
#include "pzgmesh.h"                       // for TPZGeoMesh
#include "pzgraphmesh.h"                   // for TPZGraphMesh
#include "pzlog.h"                         // for glogmutex, LOGPZ_DEBUG
#include "pzmanvector.h"                   // for TPZManVector
#include "TPZMaterial.h"                    // for TPZMaterial
#include "pzmvmesh.h"                      // for TPZMVGraphMesh
#include "pzseqsolver.h"                   // for TPZSequenceSolver
#include "pzsolve.h"                       // for TPZMatrixSolver, TPZSolver
#include "pzstack.h"                       // for TPZStack
#include "pzstepsolver.h"                  // for TPZStepSolver
#include "pzstrmatrix.h"                   // for TPZStructMatrix, TPZStruct...
#include "pzv3dmesh.h"                     // for TPZV3DGraphMesh
#include "pzvec.h"                         // for TPZVec, operator<<
#include "pzvtkmesh.h"                     // for TPZVTKGraphMesh
#include "tpznodesetcompute.h"             // for TPZNodesetCompute
#include "tpzsparseblockdiagonal.h"        // for TPZSparseBlockDiagonal

#ifdef WIN32
#include "pzsloan.h"                       // for TPZSloan
#endif

#ifdef PZ_LOG
static TPZLogger logger("pz.analysis");
static TPZLogger loggerError("pz.analysis.error");
#endif

//@orlandini: does anyone know if boost renumbering still works?
#undef USE_BOOST_RENUMBERING
#if defined(USING_BOOST) && defined(USE_BOOST_RENUMBERING)
#include "TPZBoostGraph.h"
/**
 * @brief Renumbering will use boost library.
 * @ingroup analysis
 */
#define RENUMBER TPZBoostGraph(TPZBoostGraph::KMCExpensive))
#else
/**
 * @brief Renumbering will use sloan library.
 * @ingroup analysis
 */
#define RENUMBER TPZSloanRenumbering()
//#define RENUMBER TPZCutHillMcKee()
#endif

using namespace std;

void TPZAnalysis::SetStructuralMatrix(TPZStructMatrix &strmatrix){
	fStructMatrix = TPZAutoPointer<TPZStructMatrix>(strmatrix.Clone());
}

void TPZAnalysis::SetStructuralMatrix(TPZAutoPointer<TPZStructMatrix> strmatrix){
	fStructMatrix = TPZAutoPointer<TPZStructMatrix>(strmatrix->Clone());
}
TPZAnalysis::TPZAnalysis() : TPZRegisterClassId(&TPZAnalysis::ClassId),
fGeoMesh(0), fCompMesh(0),
fRhs(false), fSolution(false),//TODOCOMPLEX:set matrix type (complex/real)
fSolver(0), fStep(0), fTime(0.), fNthreadsError(0),fStructMatrix(0), fRenumber(new RENUMBER)
, fGuiInterface(NULL), fTable() {
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
}


TPZAnalysis::TPZAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth, std::ostream &out) :
TPZRegisterClassId(&TPZAnalysis::ClassId),
fGeoMesh(0), fCompMesh(0),
fRhs(false), fSolution(false),//TODOCOMPLEX:set matrix type (complex/real)
fSolver(0), fStep(0), fTime(0.), fNthreadsError(0), fStructMatrix(0), fRenumber(new RENUMBER), fGuiInterface(NULL),  fTable()
{
	fGraphMesh[0] = 0;
	fGraphMesh[1] = 0;
	fGraphMesh[2] = 0;
    this->SetCompMesh(mesh, mustOptimizeBandwidth);
}

TPZAnalysis::TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth, std::ostream &out) :
TPZRegisterClassId(&TPZAnalysis::ClassId),
fGeoMesh(0), fCompMesh(0),
fRhs(false), fSolution(false),//TODOCOMPLEX:set matrix type (complex/real)
fSolver(0), fStep(0), fTime(0.), fNthreadsError(0),fStructMatrix(0), fRenumber(new RENUMBER), fGuiInterface(NULL),  fTable()
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
//        fCompMesh->InitializeBlock();
        int64_t neq = fCompMesh->NEquations();
        if(neq > 20000 && mustOptimizeBandwidth)
        {
            std::cout << __PRETTY_FUNCTION__ << " optimizing bandwidth\n";
            std::cout.flush();
        }
        if(mustOptimizeBandwidth)
        {
            OptimizeBandwidth();
        }
        if(neq > 20000 && mustOptimizeBandwidth)
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
    if(fSolver) {
        delete fSolver;
        fSolver = NULL;
    }
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
//    fCompMesh->InitializeBlock();
	TPZVec<int64_t> perm,iperm;
	
	TPZStack<int64_t> elgraph;
	TPZStack<int64_t> elgraphindex;
	int64_t nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int64_t nel = elgraphindex.NElements()-1;
	TPZSloan sloan(nel,nindep);
	sloan.SetElementGraph(elgraph,elgraphindex);
	sloan.Resequence(perm,iperm);
	fCompMesh->Permute(perm);
#else
	if(!fCompMesh) return;
//    fCompMesh->InitializeBlock();
	
	TPZVec<int64_t> perm,iperm;
	
	TPZStack<int64_t> elgraph,elgraphindex;
	int64_t nindep = fCompMesh->NIndependentConnects();
    
    /// if there are no connects, there is no bandwidth to be optimized
    if(nindep == 0) return;
    
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int64_t nel = elgraphindex.NElements()-1;
	int64_t el,ncel = fCompMesh->NElements();
	int maxelcon = 0;
	for(el = 0; el<ncel; el++)
	{
		TPZCompEl *cel = fCompMesh->ElementVec()[el];
		if(!cel) continue;
		std::set<int64_t> indepconlist,depconlist;
		cel->BuildConnectList(indepconlist,depconlist);
		int64_t locnindep = indepconlist.size();
		maxelcon = maxelcon < locnindep ? locnindep : maxelcon;
	}
	fRenumber->SetElementsNodes(nel,nindep);
	fRenumber->SetElementGraph(elgraph,elgraphindex);
#ifdef PZ_LOG2
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        fRenumber->Print(elgraph, elgraphindex, "Elgraph of submesh", sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    bool shouldstop = false;
    try {
        fRenumber->Resequence(perm,iperm);
    } catch(...)
    {
        fRenumber->PlotElementGroups("ElementGroups.vtk", fCompMesh);
        fCompMesh->Print();
        shouldstop = true;
    }
    if(shouldstop) DebugStop();
    
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
	int64_t sz = this->Mesh()->NEquations();
	this->Rhs().Redim(sz,numloadcases);
    //int64_t othersz = fStructMatrix->Mesh()->NEquations();
	fStructMatrix->Assemble(this->Rhs(),fGuiInterface);
}//void

template<class TVar>
void TPZAnalysis::AssembleInternal()
{
    auto mySolver =
        dynamic_cast<TPZMatrixSolver<TVar> *>(fSolver);
    if(!fCompMesh || !fStructMatrix || !mySolver)
	{
		std::stringstream sout;
		sout << "TPZAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh
		<< " fStructMatrix " << (void *) fStructMatrix.operator->()
		<< " mySolver " << (void *) mySolver;
#ifndef WINDOWS
		sout << " at file " << __FILE__ << " line " << __LINE__ ;
#else
		sout << " TPZAnalysis::Assemble() " ;
#endif
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str().c_str());
#else
		std::cout << sout.str().c_str() << std::endl;
#endif
		return;
	}
    int numloadcases = ComputeNumberofLoadCases();
	int64_t sz = fCompMesh->NEquations();
	fRhs.Redim(sz,numloadcases);
	if(mySolver->Matrix() && mySolver->Matrix()->Rows()==sz)
	{
		mySolver->Matrix()->Zero();
		fStructMatrix->Assemble(*(mySolver->Matrix().operator ->()),fRhs,fGuiInterface);
	}
	else
	{
        
		auto *mat = fStructMatrix->CreateAssemble(fRhs,fGuiInterface);
		mySolver->SetMatrix(mat);
		//aqui TPZFMatrix<TVar> nao eh nula
	}
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        PrintVectorByElement(sout, fRhs, 1.e-6);
//        fRhs.Print("Rhs",sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
	mySolver->UpdateFrom(mySolver->Matrix());
}
void TPZAnalysis::Assemble()
{
	//TODOCOMPLEX
    TPZAnalysis::AssembleInternal<STATE>();
}


template<class TVar>
void TPZAnalysis::SolveInternal(){
	int64_t numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) 
    {
        DebugStop();
    }
	int64_t nReducedEq = fStructMatrix->NReducedEquations();
    auto mySolver =
        dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
    if (nReducedEq == numeq) 
    {
        TPZFMatrix<TVar> residual(fRhs);
        TPZFMatrix<TVar> delu(numeq,1,0.);
        //      TVar normres  = Norm(residual);
        //	cout << "TPZAnalysis::Solve residual : " << normres << " neq " << numeq << endl;
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            TPZFMatrix<TVar> res2(fRhs);
            mySolver->Matrix()->Residual(fSolution,fRhs,res2);
            std::stringstream sout;
            sout << "Residual norm " << Norm(res2) << std::endl;
    //		res2.Print("Residual",sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
    
//        {
//            std::ofstream out("Matrix.nb");
//            mySolver->Matrix()->Print("Stiffness = ",out,EMathematicaInput);
//
//        }
        mySolver->Solve(residual, delu);
        fSolution = delu;
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            if(!mySolver->Matrix()->IsDecomposed())
            {
                TPZFMatrix<TVar> res2(fRhs);
                mySolver->Matrix()->Residual(delu,fRhs,res2);
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
        TPZFMatrix<TVar> residual(nReducedEq,1,0.);
    	TPZFMatrix<TVar> delu(nReducedEq,1,0.);
        fStructMatrix->EquationFilter().Gather(fRhs,residual);
	    mySolver->Solve(residual, delu);
        fSolution.Redim(numeq,1);
        fStructMatrix->EquationFilter().Scatter(delu,fSolution);
    }
#ifdef PZ_LOG
    std::stringstream sout;
    TPZStepSolver<TVar> *step = dynamic_cast<TPZStepSolver<TVar> *> (mySolver);
    if(!step) DebugStop();
    int64_t nsing = step->Singular().size();
	if(nsing && logger.isWarnEnabled()) {
		sout << "Number of singular equations " << nsing;
		std::list<int64_t>::iterator it = step->Singular().begin();
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
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
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

void TPZAnalysis::Solve() {
    SolveInternal<STATE>();
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
        int64_t i,nelements = fCompMesh->ConnectVec().NElements();
        for(i=0;i<nelements;i++) {
            TPZConnect &gnod = fCompMesh->ConnectVec()[i];
            if(gnod.SequenceNumber()!=-1) {
                out << "Connect index " << i << endl;
                gnod.Print(*fCompMesh,out);
            }
        }
    }
    TPZBaseMatrix &sol = this->Solution();
	sol.Print("fSolution",out);
    if (fCompMesh) {
        fCompMesh->ConnectSolution(out);
    }
}


void TPZAnalysis::PostProcess(TPZVec<REAL> &ervec, std::ostream &out) {
	int64_t i;
	//TPZVec<REAL> ux((int) neq);
	//TPZVec<REAL> sigx((int) neq);
	TPZManVector<REAL,10> values(10,0.);
	TPZManVector<REAL,10> values2(10,0.);
	fCompMesh->LoadSolution(fSolution);
	TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
	TPZManVector<REAL,10> errors(10);
	errors.Fill(0.0);
	int64_t nel = elvec.NElements();
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
      if (el->Material() && el->Material()->HasExactSol()) {
        el->EvaluateError(errors, 0);
      } else {
          PZError<<__PRETTY_FUNCTION__;
          PZError<<" the material has no associated exact solution\n";
          PZError<<"Aborting...";
          DebugStop();
      }
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

void TPZAnalysis::PostProcessError(TPZVec<REAL> &ervec, bool store_error, std::ostream &out ){
  if(!fNthreadsError){
    PostProcessErrorSerial(ervec, store_error, out);
  }
  else{
    PostProcessErrorParallel(ervec, store_error, out);
  }
}

void TPZAnalysis::CreateListOfCompElsToComputeError(TPZAdmChunkVector<TPZCompEl *> &elvecToComputeError){
  
  int64_t neq = fCompMesh->NEquations();
  TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
  const int64_t ncompel = elvec.NElements();
  elvecToComputeError.Resize(ncompel);
  int64_t i, nel = elvec.NElements();
  int64_t nelToCompute = 0;
  for(i=0;i<nel;i++) {
    TPZCompEl *el = (TPZCompEl *) elvec[i];
    if(el) {
      TPZMaterial *mat = el->Material();
      TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
      if(!bc){
        elvecToComputeError[nelToCompute] = el;
        nelToCompute++;
      }
    }//if(el)
  }//i
  
  elvecToComputeError.Resize(nelToCompute);
  
}

void *TPZAnalysis::ThreadData::ThreadWork(void *datavoid)
{
  ThreadData *data = (ThreadData *) datavoid;
  const int64_t nelem = data->fElvec.NElements();
  TPZManVector<REAL,10> errors(10);
 
  // Getting unique id for each thread
  const int64_t myid = [&]()
  {
    std::scoped_lock lock(data->fMutexThreadId);
    const int64_t myid_loc = data->ftid;
    data->ftid++;
    return myid_loc;
  }();

  
  
  do{

    const int64_t iel = [&]()
    {
      std::scoped_lock lock(data->fMutexAccessEl);
      const int64_t iel_loc = data->fNextElement;
      data->fNextElement++;
      return iel_loc;
    }();
    
    // For all the elements it tries to get after the last one
    if ( iel >= nelem ) continue;
    
    TPZCompEl *cel = data->fElvec[iel];
    if (cel->Material() && cel->Material()->HasExactSol()) {
      cel->EvaluateError(errors, data->fStoreError);
    } else {
      PZError<<__PRETTY_FUNCTION__;
      PZError<<" the material has no associated exact solution\n";
      PZError<<"Aborting...";
      DebugStop();
    }
    
    const int nerrors = errors.NElements();
    data->fvalues[myid].Resize(nerrors, 0.);


    for(int ier = 0; ier < nerrors; ier++)
    {
      (data->fvalues[myid])[ier] += errors[ier] * errors[ier];
    }
    
    
  } while (data->fNextElement < nelem);
  return data;
}

void TPZAnalysis::PostProcessErrorParallel(TPZVec<REAL> &ervec, bool store_error, std::ostream &out ){ //totto

  fCompMesh->LoadSolution(fSolution);
  const int numthreads = this->fNthreadsError;
  std::vector<std::thread> allthreads;

  TPZAdmChunkVector<TPZCompEl *> elvec;
  std::chrono::high_resolution_clock::now();
  auto tsim1 = std::chrono::high_resolution_clock::now();
  CreateListOfCompElsToComputeError(elvec);
  auto tsim2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, milli> elapsed_time = tsim2 - tsim1;
  std::cout << "Total wall time of CreateListOfCompElsToComputeError = " << elapsed_time.count()<<"ms"<<std::endl;
  
  
  ThreadData threaddata(elvec,store_error);
  threaddata.fvalues.Resize(numthreads);
  for(int iv = 0 ; iv < numthreads ; iv++){
      threaddata.fvalues[iv].Resize(10);
      threaddata.fvalues[iv].Fill(0.0);
  }


  auto tthread1 = std::chrono::high_resolution_clock::now();
  
  for(int itr=0; itr<numthreads; itr++)
  {
    allthreads.push_back(std::thread(ThreadData::ThreadWork, &threaddata));
  }
  
  for(int itr=0; itr<numthreads; itr++)
  {
    allthreads[itr].join();
  }

  auto tthread2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, milli> elapsed_time_thread = tthread2 - tthread1;
  std::cout << "Total wall time of ThreadWork = " << elapsed_time_thread.count()<<"ms"<<std::endl;
  
  
  // Sanity check. There should be number of ids equal to number of threads
  if(threaddata.ftid != numthreads){
    DebugStop();
  }

  
  TPZManVector<REAL,10> values;
  // Assuming the first is equal to the others
  const int nerrors = threaddata.fvalues[0].NElements();
  values.Resize(nerrors,0);
  // Summing up all the values of all threads
  for(int it = 0 ; it < numthreads ; it++){
      for(int ir = 0 ; ir < nerrors ; ir++){
          values[ir] += (threaddata.fvalues[it])[ir];
      }
  }
  
  //const int nerrors = values.NElements();
  ervec.Resize(nerrors);
  ervec.Fill(-10.0);
    
  if (nerrors < 3) {
    PZError << endl << "TPZAnalysis::PostProcess - At least 3 norms are expected." << endl;
    out<<endl<<"############"<<endl;

#ifdef PZDEBUG
    for(int ier = 0; ier < nerrors; ier++)
      out << endl << "error " << ier << "  = " << sqrt(values[ier]);
#endif
  }
  else{
#ifdef PZDEBUG
    out << "############" << endl;
    out <<"Norma H1 or L2 -> p = "  << sqrt(values[0]) << endl;
    out <<"Norma L2 or L2 -> u = "    << sqrt(values[1]) << endl;
    out << "Semi-norma H1 or L2 -> div = "    << sqrt(values[2])  <<endl;
    for(int ier = 3; ier < nerrors; ier++)
      out << "other norms = " << sqrt(values[ier]) << endl;
#endif
  }
  
  // Returns the calculated errors.
  for(int i=0;i<nerrors;i++)
    ervec[i] = sqrt(values[i]);
  return;
  
}

#include "pzsubcmesh.h"

void TPZAnalysis::PostProcessErrorSerial(TPZVec<REAL> &ervec, bool store_error, std::ostream &out ){

    fCompMesh->EvaluateError(store_error, ervec);


    const int nerrors = ervec.NElements();

    if (nerrors < 3) {
        PZError << endl << "TPZAnalysis::PostProcess - At least 3 norms are expected." << endl;
        out<<endl<<"############"<<endl;
        for(int ier = 0; ier < nerrors; ier++)
            out << endl << "error " << ier << "  = " << ervec[ier];
    }
    else {
#ifdef PZDEBUG
        out << "############" << endl;
        out << "Norma H1 or L2 -> p = " << ervec[0] << std::endl;
        out << "Norma L2 or L2 -> u = " << ervec[1] << std::endl;
        out << "Semi-norma H1 or L2 -> div = " << ervec[2] << std::endl;
        for (int ier = 3; ier < ervec.size(); ier++)
            out << "other norms = " << ervec[ier] << std::endl;
#endif
    }
}

void TPZAnalysis::PostProcessTable( TPZFMatrix<REAL> &,std::ostream & )//pos,out
{
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	int64_t nel = elvec.NElements();
	for(int64_t i=0;i<nel;i++) {
		//TPZCompEl *el = (TPZCompEl *) elvec[i];
		//if(el)	el->PrintTable(fExact,pos,out);
	}
	return;
}

void TPZAnalysis::SetExact(std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)>f,int p){
    //TODOCOMPLEX
    TPZAnalysis::SetExactInternal<STATE>(f,p);
}

template<class TVar>
void TPZAnalysis::SetExactInternal(std::function<void (const TPZVec<REAL> &loc, TPZVec<TVar> &result, TPZFMatrix<TVar> &deriv)>f,int p){
    if(!fCompMesh){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" there is no associated computational mesh.\n";
        PZError<<"Aborting...\n";
        DebugStop();
    }
    for(auto imat : fCompMesh->MaterialVec()){
        imat.second->SetExactSol(f,p);
    }
}
template<class TVar>
void TPZAnalysis::ShowShapeInternal(
    const TPZStack<std::string> &scalnames,
    const TPZStack<std::string> &vecnames,
    const std::string &plotfile,
    TPZVec<int64_t> &equationindices){
  DefineGraphMesh(fCompMesh->Dimension(), scalnames, vecnames, plotfile);
  int porder = fCompMesh->GetDefaultOrder();

  int neq = equationindices.size();
  TPZFMatrix<TVar> solkeep(fSolution);
  TPZFMatrix<TVar> &mysol = fSolution;
  mysol.Zero();
  for (int ieq = 0; ieq < neq; ieq++) {
    mysol(equationindices[ieq], 0) = 1.;
    LoadSolution();
    Mesh()->TransferMultiphysicsSolution();
    PostProcess(porder + 1);
    fSolution.Zero();
  }
  fSolution = solkeep;
  LoadSolution();
}
void TPZAnalysis::ShowShape(const std::string &plotfile, TPZVec<int64_t> &equationindices)
{
	
    SetStep(1);
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = fCompMesh->MaterialVec().begin()->second;
    int nstate = mat->NStateVariables();
    if (nstate == 1) {
        scalnames.Push("Solution");
    }
    else
    {
        vecnames.Push("Solution");
    }
    //TODOCOMPLEX
    ShowShapeInternal<STATE>(scalnames,vecnames,plotfile,equationindices);
}

void TPZAnalysis::ShowShape(const std::string &plotfile, TPZVec<int64_t> &equationindices, int matid, const TPZVec<std::string> &varname)
{
    TPZMaterial *mat = fCompMesh->FindMaterial(matid);
    if(!mat) DebugStop();
    
    int n_varnames = varname.size();
    TPZStack<std::string> scalnames,vecnames;
    for(int i_var = 0 ; i_var<n_varnames ; i_var++){
        int varindex = mat->VariableIndex(varname[i_var]);
        if(varindex == -1) DebugStop();
        SetStep(1);
        
        int nstate = mat->NSolutionVariables(varindex);
        if (nstate == 1) {
            scalnames.Push(varname[i_var]);
        }
        else
        {
            vecnames.Push(varname[i_var]);
        }
    }
    //TODOCOMPLEX
    ShowShapeInternal<STATE>(scalnames,vecnames,plotfile,equationindices);
}



void TPZAnalysis::LoadShape(double ,double , int64_t ,TPZConnect* start){
	//void TPZAnalysis::LoadShape(double dx,double dy, int numelem,TPZConnect* start){
	Assemble();
	fRhs.Zero();
	fSolution.Zero();
	
}


void TPZAnalysis::Run(std::ostream &out)
{
    int64_t neq = fCompMesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
	Assemble();

    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    Solve();
    const auto t3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, milli> time_assemble = t2-t1;
    std::chrono::duration<double, milli> time_solve = t3-t2;
    std::cout << "Time for assembly " << time_assemble.count() << "ms Time for solving " << time_solve.count() << "ms"<<std::endl;
}

/** @brief Define GrapMesh as V3D, DX, MV or VTK depending on extension of the file */
void TPZAnalysis::DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile){
    std::set<int> matids;
    IdentifyPostProcessingMatIds(dimension, matids);
    this->DefineGraphMesh(dimension, matids, scalnames, vecnames, plotfile);
    
}
/** @brief Define GrapMesh as VTK with tensorial names depending on extension of the file */
void TPZAnalysis::DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile){
    std::set<int> matids;
    IdentifyPostProcessingMatIds(dimension, matids);
    this->DefineGraphMesh(dimension, matids, scalnames, vecnames, tensnames, plotfile);
}

void TPZAnalysis::DefineGraphMesh(int dim, const std::set<int> & matids , const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile) {

    int dim1 = dim - 1;
    if (!fCompMesh || matids.size() == 0) {
        cout << "TPZAnalysis::DefineGraphMesh nothing to post-process\n";
        return;
    }

    if (fGraphMesh[dim1]) delete fGraphMesh[dim1];
    fScalarNames[dim1] = scalnames;
    fVectorNames[dim1] = vecnames;
    int posplot = plotfile.rfind(".plt");
    int64_t filelength = plotfile.size();
    if (filelength - posplot == 3) {
        if(tensnames.size() != 0)
        {
            std::cout << "Tensors are not beint post-process for file = " << plotfile << std::endl;
        }
        fGraphMesh[dim1] = new TPZV3DGraphMesh(fCompMesh, dim, matids, scalnames, vecnames);
    } else {
        int posdx = plotfile.rfind(".dx");
        if (filelength - posdx == 3) {
            if(tensnames.size() != 0) {
                std::cout << "Tensors are not beint post-process for file = " << plotfile << std::endl;
                
            }
            fGraphMesh[dim1] = new TPZDXGraphMesh(fCompMesh, dim, matids, scalnames, vecnames);
        } else {
            int pospos = plotfile.rfind(".pos");
            if (filelength - pospos == 3) {
                if(tensnames.size() != 0) {
                    std::cout << "Tensors are not beint post-process for file = " << plotfile << std::endl;
                }
                fGraphMesh[dim1] = new TPZMVGraphMesh(fCompMesh, dim, matids, scalnames, vecnames);
            } else {
                int posvtk = plotfile.rfind(".vtk");
                if (filelength - posvtk == 4) {
                    fGraphMesh[dim1] = new TPZVTKGraphMesh(fCompMesh, dim, matids, scalnames, vecnames, tensnames);
                } else {
                    cout << "grafgrid was not created\n";
                    fGraphMesh[dim1] = 0;
                }
            }
        }
    }
    if (fGraphMesh[dim1]) {
        fGraphMesh[dim1]->SetFileName(plotfile);
    }
}

void TPZAnalysis::DefineGraphMesh(int dim, const std::set<int> & matids, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile){
    TPZVec<std::string> tensnames;
    this->DefineGraphMesh(dim,matids, scalnames, vecnames, tensnames, plotfile);
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

/** @brief Fill mat ids with materials with dimension dim wich are not boundary conditions or interface  */
void TPZAnalysis::IdentifyPostProcessingMatIds(int dimension, std::set<int> & matids){
    std::map<int, TPZMaterial * >::iterator matit;
    for (matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++) {
        
        TPZMaterial *mat = matit->second;
        TPZBndCond *bc = dynamic_cast<TPZBndCond *> (mat);
        TPZLagrangeMultiplier *lag = dynamic_cast<TPZLagrangeMultiplier *> (mat);
        if (mat && !bc && !lag && mat->Dimension() == dimension){
            matids.insert(mat->Id());
        }
    }
    if (matids.size() == 0) {
        std::cout << __PRETTY_FUNCTION__ << " No materials were identified to perform post-processing. \n";
        DebugStop();
        return;
    }
}

void TPZAnalysis::PostProcess(int resolution, int dimension){
	int dim1 = dimension-1;
	if(!fGraphMesh[dim1]) return;

	fGraphMesh[dim1]->SetCompMesh(fCompMesh,fGraphMesh[dim1]->MaterialIds());
	fGraphMesh[dim1]->SetResolution(resolution);
	fGraphMesh[dim1]->DrawMesh(1);
	fGraphMesh[dim1]->DrawSolution(fStep,fTime);
	fStep++;
}


template<class TVar>
void TPZAnalysis::AnimateRunInternal(
    int64_t num_iter, int steps, TPZVec<std::string> &scalnames,
    TPZVec<std::string> &vecnames, const std::string &plotfile)
{
    Assemble();
	int64_t numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return;
	
	TPZFMatrix<TVar> residual(fRhs);
	int dim = HighestDimension();
    std::set<int> matids;
    IdentifyPostProcessingMatIds(dim, matids);
	TPZDXGraphMesh gg(fCompMesh,dim,matids,scalnames,vecnames) ;
	gg.SetFileName(plotfile);
	gg.SetResolution(0);
    auto mySolver = dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
	for(auto i=1; i<=num_iter;i+=steps){
		TPZStepSolver<TVar> sol;
		sol.ShareMatrix(*mySolver);
		sol.SetJacobi(i,0.,0);
		SetSolver(sol);
		mySolver->Solve(fRhs, fSolution);
		fCompMesh->LoadSolution(fSolution);
		gg.DrawSolution(i-1,0);
	}
}
    
void TPZAnalysis::AnimateRun(
    int64_t num_iter, int steps, TPZVec<std::string> &scalnames,
    TPZVec<std::string> &vecnames, const std::string &plotfile)
{
	//TODOCOMPLEX
    AnimateRunInternal<STATE>(num_iter,steps,scalnames,vecnames,plotfile);
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
fGeoElId(), fCompElPtr(), fDimension(-1), fLocations(), fVariableNames() {
}

int TPZAnalysis::TTablePostProcess::ClassId() const{
    return Hash("TPZAnalysis::TTablePostProcess");
}

void TPZAnalysis::TTablePostProcess::Write(TPZStream &buf, int withclassid) const {
    buf.Write(fGeoElId);
    buf.WritePointers(fCompElPtr);
    buf.Write(&fDimension);
    buf.Write(fLocations);
    buf.Write(fVariableNames);
}

void TPZAnalysis::TTablePostProcess::Read(TPZStream &buf, void *context){
    buf.Read(fGeoElId);
    buf.ReadPointers(fCompElPtr);
    buf.Read(&fDimension);
    buf.Read(fLocations);
    buf.Read(fVariableNames);
}


TPZAnalysis::TTablePostProcess::~TTablePostProcess() {
	fDimension = -1;
}

void TPZAnalysis::DefineElementTable(int dimension, TPZVec<int64_t> &GeoElIds, TPZVec<REAL> &points) {
	fTable.fDimension = dimension;
	fTable.fGeoElId = GeoElIds;
	fTable.fLocations = points;
	int64_t numel = GeoElIds.NElements();
	fTable.fCompElPtr.Resize(numel);
	int64_t iel;
	for(iel=0; iel<numel; iel++) {
		TPZGeoEl *gel = (TPZGeoEl *) fGeoMesh->FindElement(GeoElIds[iel]);
		fTable.fCompElPtr[iel] = (gel) ? gel->Reference() : 0;
	}
}

void TPZAnalysis::SetTableVariableNames(TPZVec<std::string> varnames) {
    fTable.fVariableNames = varnames;
}

void TPZAnalysis::PrePostProcessTable(std::ostream &out_file){
	TPZCompEl *cel;
	int numvar = fTable.fVariableNames.NElements();
	for(int iv=0; iv<numvar; iv++) {
		int64_t numel = fTable.fCompElPtr.NElements();
		for(int64_t iel=0; iel<numel; iel++) {
			cel = (TPZCompEl *) fTable.fCompElPtr[iel];
			if(cel) cel->PrintTitle(fTable.fVariableNames[iv].c_str(),out_file);
		}
	}
	out_file << endl;
	int dim;
	TPZVec<REAL> point(fTable.fDimension);
	for(dim=1; dim<fTable.fDimension+1; dim++) {
		for(int iv=0; iv<numvar; iv++) {
			int64_t numel = fTable.fCompElPtr.NElements();
			for(int64_t iel=0; iel<numel; iel++) {
				int d;
				for(d=0; d<fTable.fDimension; d++) {
					point[d] = fTable.fLocations[iel*fTable.fDimension+d];
				}
				cel = (TPZCompEl *) fTable.fCompElPtr[iel];
				if(cel) cel->PrintCoordinate(point,dim,out_file);
			}
		}
		out_file << endl;
	}
}

void TPZAnalysis::PostProcessTable(std::ostream &out_file) {
	TPZVec<REAL> point(fTable.fDimension);
	int numvar = fTable.fVariableNames.NElements();
	TPZCompEl *cel;
	for(int iv=0; iv<numvar; iv++) {
		int64_t numel = fTable.fCompElPtr.NElements();
		for(int64_t iel=0; iel<numel; iel++) {
			int d;
			for(d=0; d<fTable.fDimension; d++) {
				point[d] = fTable.fLocations[iel*fTable.fDimension+d];
			}
			cel = (TPZCompEl *) fTable.fCompElPtr[iel];
			if(cel) cel->PrintSolution(point,fTable.fVariableNames[iv].c_str(),out_file);
		}
	}
	out_file << endl;
}
void TPZAnalysis::SetSolver(const TPZSolver &solver){
	if(fSolver) delete fSolver;
    //TODOCOMPLEX
    fSolver = (TPZMatrixSolver<STATE> *) solver.Clone();
}

template<class TVar>
TPZMatrixSolver<TVar> *TPZAnalysis::BuildPreconditioner(EPrecond preconditioner, bool overlap)
{
    auto mySolver = dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
	if(!mySolver || !mySolver->Matrix())
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
		TPZStack<int64_t> elementgraph,elementgraphindex;
		int64_t nindep = fCompMesh->NIndependentConnects();
		int64_t neq = fCompMesh->NEquations();
		fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
		int64_t nel = elementgraphindex.NElements()-1;
		TPZRenumbering renum(nel,nindep);
		renum.ConvertGraph(elementgraph,elementgraphindex,nodeset.Nodegraph(),nodeset.Nodegraphindex());
		nodeset.AnalyseGraph();

		TPZStack<int64_t> blockgraph,blockgraphindex;
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
		TPZStack<int64_t> expblockgraph,expblockgraphindex;
		
		nodeset.ExpandGraph(blockgraph,blockgraphindex,fCompMesh->Block(),expblockgraph,expblockgraphindex);
#ifdef PZ_LOG
#ifdef PZDEBUG2
        if (logger.isDebugEnabled())
        {
            std::map<int64_t,int64_t> blocksizes;
            int64_t i;
            int64_t totalsize = 0;
            for(i=0; i< expblockgraphindex.NElements()-1;i++)
            {
                int64_t bls = expblockgraphindex[i+1]-expblockgraphindex[i];
                blocksizes[bls]++;
                totalsize += bls*bls;
            }
            std::map<int64_t,int64_t>::iterator it;
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
			TPZSparseBlockDiagonal<TVar> *sp = new TPZSparseBlockDiagonal<TVar>(expblockgraph,expblockgraphindex,neq);
			TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(sp);
			step->SetDirect(ELU);
			step->SetReferenceMatrix(mySolver->Matrix());
			return step;
		}
		else if (overlap)
		{
			TPZBlockDiagonalStructMatrix blstr(fCompMesh);
			TPZBlockDiagonal<TVar> *sp = new TPZBlockDiagonal<TVar>();
			blstr.AssembleBlockDiagonal(*sp);
			TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(sp);
			step->SetDirect(ELU);
			return step;
		}
		else
		{
			TPZVec<int> blockcolor;
			int numcolors = nodeset.ColorGraph(expblockgraph,expblockgraphindex,neq,blockcolor);
			return BuildSequenceSolver<TVar>(expblockgraph,expblockgraphindex,neq,numcolors,blockcolor);
		}
	}
	return 0;
}

/** @brief Build a sequence solver based on the block graph and its colors */
template<class TVar>
TPZMatrixSolver<TVar> *TPZAnalysis::BuildSequenceSolver(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex, int64_t neq, int numcolors, TPZVec<int> &colors)
{
	TPZVec<TPZMatrix<TVar> *> blmat(numcolors);
	TPZVec<TPZStepSolver<TVar> *> steps(numcolors);
	int c;
    auto mySolver =
        dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
	for(c=0; c<numcolors; c++)
	{
		blmat[c] = new TPZSparseBlockDiagonal<TVar>(graph,graphindex, neq, c, colors);
		steps[c] = new TPZStepSolver<TVar>(blmat[c]);
		steps[c]->SetDirect(ELU);
		steps[c]->SetReferenceMatrix(mySolver->Matrix());
	}
	if(numcolors == 1) return steps[0];
	TPZSequenceSolver<TVar> *result = new TPZSequenceSolver<TVar>;
	result->ShareMatrix(*mySolver);
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

template<class TVar>
TPZMatrixSolver<TVar> &TPZAnalysis::MatrixSolver(){
    const auto tmp = dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
    if(fSolver && !tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible Solver type! Aborting\n";
        DebugStop();
    }
    return *tmp;
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
    int64_t nelem = fCompMesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
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
static void ConnectSolution(int64_t cindex, TPZCompMesh *cmesh, TPZFMatrix<STATE> &glob, TPZVec<STATE> &sol)
{
    int64_t seqnum = cmesh->ConnectVec()[cindex].SequenceNumber();
    int blsize = cmesh->Block().Size(seqnum);
    int position = cmesh->Block().Position(seqnum);
    sol.resize(blsize);
    for (int64_t i=position; i< position+blsize; i++) {
        sol[i-position] = glob(i,0);
    }
}

static STATE ConnectNorm(int64_t cindex, TPZCompMesh *cmesh, TPZFMatrix<STATE> &glob)
{
    TPZManVector<STATE,20> cvec;
    ConnectSolution(cindex, cmesh, glob, cvec);
    STATE norm = 0.;
    for (int64_t i=0; i<cvec.size(); i++) {
        norm += cvec[i]*cvec[i];
    }
    norm = sqrt(norm);
    return norm;
}

/// Print the residual vector for those elements with entry above a given tolerance
void TPZAnalysis::PrintVectorByElement(std::ostream &out, TPZFMatrix<STATE> &vec, REAL tol)
{
    int64_t nel = fCompMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCompMesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        int ic;
        for (ic=0; ic<nc; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
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
            int64_t cindex = cel->ConnectIndex(ic);
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

int TPZAnalysis::ClassId() const{
    return Hash("TPZAnalysis");
}

void TPZAnalysis::Write(TPZStream &buf, int withclassid) const{
    TPZPersistenceManager::WritePointer(fGeoMesh, &buf);
    TPZPersistenceManager::WritePointer(fCompMesh, &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh[0], &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh[1], &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh[2], &buf);
    fRhs.Write(buf,withclassid);
    fSolution.Write(buf,withclassid);
    TPZPersistenceManager::WritePointer(fSolver, &buf);
    buf.Write(fScalarNames[0]);
    buf.Write(fScalarNames[1]);
    buf.Write(fScalarNames[2]);
    buf.Write(fVectorNames[0]);
    buf.Write(fVectorNames[1]);
    buf.Write(fVectorNames[2]);
    buf.Write(&fStep);
    buf.Write(&fTime);
    buf.Write(&fNthreadsError);
    TPZPersistenceManager::WritePointer(fStructMatrix.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fRenumber.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fGuiInterface.operator ->(), &buf);
    fTable.Write(buf,withclassid);
    buf.Write(fTensorNames[0]);
    buf.Write(fTensorNames[1]);
    buf.Write(fTensorNames[2]);
}

void TPZAnalysis::Read(TPZStream &buf, void *context){
    fGeoMesh = dynamic_cast<TPZGeoMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fCompMesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh[0] = dynamic_cast<TPZGraphMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh[1] = dynamic_cast<TPZGraphMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh[2] = dynamic_cast<TPZGraphMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fRhs.Read(buf,context);
    fSolution.Read(buf,context);
    fSolver = dynamic_cast<TPZMatrixSolver<STATE>*>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(fScalarNames[0]);
    buf.Read(fScalarNames[1]);
    buf.Read(fScalarNames[2]);
    buf.Read(fVectorNames[0]);
    buf.Read(fVectorNames[1]);
    buf.Read(fVectorNames[2]);
    buf.Read(&fStep);
    buf.Read(&fTime);
    buf.Read(&fNthreadsError);
    fStructMatrix = TPZAutoPointerDynamicCast<TPZStructMatrix>(TPZPersistenceManager::GetAutoPointer(&buf));
    fRenumber = TPZAutoPointerDynamicCast<TPZRenumbering>(TPZPersistenceManager::GetAutoPointer(&buf));
    fGuiInterface = TPZAutoPointerDynamicCast<TPZGuiInterface>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTable.Read(buf,context);
    buf.Read(fTensorNames[0]);
    buf.Read(fTensorNames[1]);
    buf.Read(fTensorNames[2]);
}

TPZAnalysis::ThreadData::ThreadData(TPZAdmChunkVector<TPZCompEl *> &elvec, bool store_error) : fNextElement(0), fvalues(0), fStoreError(store_error), ftid(0), fElvec(elvec){
}


TPZAnalysis::ThreadData::~ThreadData()
{
}

template
TPZMatrixSolver<STATE> *TPZAnalysis::BuildPreconditioner<STATE>(
    EPrecond preconditioner,bool overlap);
template
TPZMatrixSolver<STATE> &TPZAnalysis::MatrixSolver<STATE>();
