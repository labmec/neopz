#include "TPZLinearAnalysis.h"
#include "pzcmesh.h"
#include "Hash/TPZHash.h"
#include "TPZMatrixSolver.h"
#include "TPZSkylineNSymStructMatrix.h"    // for TPZSkylineNSymStructMatrix
#include "pzstepsolver.h"
#include "pzdxmesh.h"                      // for TPZDXGraphMesh
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.analysis");
static TPZLogger loggerError("pz.analysis.error");
#endif

TPZLinearAnalysis::TPZLinearAnalysis() : TPZAnalysis()
{
}

TPZLinearAnalysis::TPZLinearAnalysis(TPZCompMesh *mesh,
                                     bool mustOptimizeBandwidth,
                                     std::ostream &out) :
  TPZAnalysis(mesh,mustOptimizeBandwidth,out),
  fRhs(fSolType == EComplex ? true : false)
{
}

TPZLinearAnalysis::TPZLinearAnalysis(TPZAutoPointer<TPZCompMesh> mesh,
                                     bool mustOptimizeBandwidth,
                                     std::ostream &out) :
  TPZAnalysis(mesh,mustOptimizeBandwidth,out),
  fRhs(fSolType == EComplex ? true : false)
{
}

void TPZLinearAnalysis::Assemble()
{
  if(fSolType == EReal)
    AssembleT<STATE>();
  else
    AssembleT<CSTATE>();
}

template<class TVar>
void TPZLinearAnalysis::AssembleT()
{
  if(!fCompMesh){
    std::stringstream sout;
    sout<<__PRETTY_FUNCTION__;
    sout<<"\nNo computational mesh found!\n";
#ifdef PZ_LOG
    LOGPZ_ERROR(logger,sout.str().c_str());
#else
    std::cout << sout.str().c_str() << std::endl;
#endif
    return;
  }
  if(!this->fStructMatrix){
    std::cout<<"Setting default struct matrix: skyline(non-symmetric)"<<std::endl;
    TPZSkylineNSymStructMatrix<TVar> defaultMatrix(fCompMesh);
    this->SetStructuralMatrix(defaultMatrix);
  }
  if(!fSolver){
    std::cout<<"Setting default solver: LU"<<std::endl;
    TPZStepSolver<TVar> defaultSolver;
    defaultSolver.SetDirect(ELU);
    this->SetSolver(defaultSolver);
  }
  auto mySolver =
    dynamic_cast<TPZMatrixSolver<TVar> *>(fSolver);
  
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

void TPZLinearAnalysis::AssembleResidual()
{
  int numloadcases = ComputeNumberofLoadCases();

  if(Mesh()->NEquations()==0)
    {
      PZError << __PRETTY_FUNCTION__;
      PZError << "\nERROR: Null mesh.\nAborting...\n";
      DebugStop();
    }
    if (!this->fStructMatrix) {
      if(fSolType == EReal){
        TPZSkylineNSymStructMatrix<STATE> defaultMatrix(fCompMesh);
        this->SetStructuralMatrix(defaultMatrix);
      }
      else{
        TPZSkylineNSymStructMatrix<CSTATE> defaultMatrix(fCompMesh);
        this->SetStructuralMatrix(defaultMatrix);
      }
    }
    int64_t sz = this->Mesh()->NEquations();
	this->Rhs().Redim(sz,numloadcases);
  //int64_t othersz = fStructMatrix->Mesh()->NEquations();
	fStructMatrix->Assemble(this->Rhs(),fGuiInterface);
}//void

void TPZLinearAnalysis::Solve()
{
  if(fSolType == EReal)
    SolveT<STATE>();
  else
    SolveT<CSTATE>();
}

template<class TVar>
void TPZLinearAnalysis::SolveT()
{
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
      REAL resnorm = Norm(residual);
      if(IsZero(resnorm))
        {
          delu.Zero();
        }
      else
        {
          mySolver->Solve(residual, delu);
        }
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

void TPZLinearAnalysis::AnimateRun(
    int64_t num_iter, int steps, TPZVec<std::string> &scalnames,
    TPZVec<std::string> &vecnames, const std::string &plotfile)
{
	//TODOCOMPLEX
  AnimateRunT<STATE>(num_iter,steps,scalnames,vecnames,plotfile);
}

template<class TVar>
void TPZLinearAnalysis::AnimateRunT(
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

void TPZLinearAnalysis::SetSolver(const TPZSolver &solver){
	if(fSolver) delete fSolver;
    auto *tmpState =
      dynamic_cast<const TPZMatrixSolver<STATE>*>(&solver);
    auto *tmpCState =
      dynamic_cast<const TPZMatrixSolver<CSTATE>*>(&solver);
    if(tmpState && fSolType == EReal){
      fSolver = (TPZMatrixSolver<STATE> *) solver.Clone();
      return;
    }
    else if(tmpCState && fSolType == EComplex){
      fSolver = (TPZMatrixSolver<CSTATE> *) solver.Clone();
      return;
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" Incompatible types!\n";
    PZError<<" Aborting...\n";
    DebugStop();
}


template<class TVar>
TPZMatrixSolver<TVar> &TPZLinearAnalysis::MatrixSolver(){
    const auto tmp = dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
    if(fSolver && !tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible Solver type! Aborting\n";
        DebugStop();
    }
    return *tmp;
}

void TPZLinearAnalysis::PostProcess(int resolution, int dimension){
	int dim1 = dimension-1;
	if(!fGraphMesh[dim1]) return;

	fGraphMesh[dim1]->SetCompMesh(fCompMesh,fGraphMesh[dim1]->MaterialIds());
	fGraphMesh[dim1]->SetResolution(resolution);
	fGraphMesh[dim1]->DrawMesh(1);
	fGraphMesh[dim1]->DrawSolution(fStep,fTime);
	fStep++;
}


int TPZLinearAnalysis::ClassId() const
{
  return Hash("TPZLinearAnalysis") ^
    TPZAnalysis::ClassId() << 1;
}
  
void TPZLinearAnalysis::Write(TPZStream &buf, int withclassid) const
{
  TPZAnalysis::Write(buf,withclassid);
  fRhs.Write(buf,withclassid);
  buf.Write(&fTime);
}

void TPZLinearAnalysis::Read(TPZStream &buf, void *context)
{
  TPZAnalysis::Read(buf,context);
  fRhs.Read(buf,context);
  buf.Read(&fTime);
}

#define INSTANTIATE_TEMPLATES(TVar)                                 \
  template                                                          \
  TPZMatrixSolver<TVar> &TPZLinearAnalysis::MatrixSolver<TVar>();



INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)


#undef INSTANTIATE_TEMPLATES

template class TPZRestoreClass<TPZLinearAnalysis>;