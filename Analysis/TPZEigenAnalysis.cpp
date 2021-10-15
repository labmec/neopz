//
//
// Created by Francisco Teixeira Orlandini on 11/17/17.

#include "TPZEigenAnalysis.h"
#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZSpStructMatrix.h"
#include "pzysmp.h"
#include "pzsysmp.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZMatGeneralisedEigenVal.h"
#include "TPZMaterial.h"
#include "pzcmesh.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.analysis");
static TPZLogger loggerError("pz.analysis.error");
#endif

TPZEigenAnalysis::TPZEigenAnalysis()
  : TPZRegisterClassId(&TPZEigenAnalysis::ClassId),TPZAnalysis()
{

}

TPZEigenAnalysis::TPZEigenAnalysis(TPZCompMesh *mesh,
                                   bool mustOptimizeBandwidth, std::ostream &out)
  : TPZAnalysis(mesh, mustOptimizeBandwidth,out)
{

}

TPZEigenAnalysis::TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh,
                                   bool mustOptimizeBandwidth,
                                   std::ostream &out)
  : TPZAnalysis(mesh, mustOptimizeBandwidth, out)
{

}

template<class TVar>
TPZEigenSolver<TVar> &TPZEigenAnalysis::EigenSolver()
{
    const auto tmp = dynamic_cast<TPZEigenSolver<TVar>*>(fSolver);
    if(fSolver && !tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible Solver type! Aborting\n";
        DebugStop();
    }
    return *tmp;
}

void TPZEigenAnalysis::SetSolver(const TPZSolver &solver)
{
    if(fSolver) delete fSolver;
    auto *tmpState =
      dynamic_cast<const TPZEigenSolver<STATE>*>(&solver);
    auto *tmpCState =
      dynamic_cast<const TPZEigenSolver<CSTATE>*>(&solver);
    if(tmpState && fSolType == EReal){
      fSolver = (TPZEigenSolver<STATE> *) solver.Clone();
      return;
    }
    else if(tmpCState && fSolType == EComplex){
      fSolver = (TPZEigenSolver<CSTATE> *) solver.Clone();
      return;
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" Incompatible types!\n";
    PZError<<" Aborting...\n";
    DebugStop();
}

void TPZEigenAnalysis::Assemble()
{
  if(fSolType == EReal)
    AssembleT<STATE>();
  else
    AssembleT<CSTATE>();
}

template<class TVar>
void TPZEigenAnalysis::AssembleT()
{
  //this lambda will be used to prevent pardiso for disturbing pivots
  auto PardisoConfig = [](TPZAutoPointer<TPZMatrix<TVar>> mat){
    auto ssp = TPZAutoPointerDynamicCast<TPZSYsmpMatrix<TVar>>(mat);
    auto sp = TPZAutoPointerDynamicCast<TPZFYsmpMatrix<TVar>>(mat);
    TPZPardisoSolver<TVar> *pardiso{nullptr};
    if(ssp){
      pardiso = &(ssp->GetPardisoControl());
    }else if(sp){
      pardiso = &(sp->GetPardisoControl());
    }
    if(pardiso){
      auto param = pardiso->GetParam();
      //do not perturb pivot elements
      param[0] = 0;
    }
  };
  
  //it wont be resized or anything.
  TPZFMatrix<TVar> dummyRhs;
  if (!fCompMesh) {
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__;
    sout << "\nNo computational mesh found!\n";
#ifdef PZ_LOG
    LOGPZ_ERROR(logger, sout.str().c_str());
#else
    std::cout << sout.str().c_str() << std::endl;
#endif
    return;
  }
  if (!this->fStructMatrix) {
#ifdef USING_MKL
    std::cout << "Setting default struct matrix: sparse(non-symmetric)"
              << std::endl;
    TPZSpStructMatrix<TVar> defaultMatrix(fCompMesh);
#else
    std::cout << "Setting default struct matrix: skyline(non-symmetric)"
              << std::endl;
    TPZSkylineNSymStructMatrix<TVar> defaultMatrix(fCompMesh);
#endif
    this->SetStructuralMatrix(defaultMatrix);
  }
  
  fStructMatrix->SetComputeRhs(false);
  
  {
    auto *eigSolver = &(this->EigenSolver<TVar>()); 
    if (!eigSolver) {
      std::cout << "Setting default solver: Krylov" << std::endl;
      constexpr int nev{10};
      constexpr int dimKrylov{100};
      TPZKrylovEigenSolver<TVar> defaultSolver;
      defaultSolver.SetNEigenpairs(nev);
      defaultSolver.SetKrylovDim(dimKrylov);
      std::cout << "Setting nev: " << nev << std::endl;
      std::cout << "Setting krylov dim: " << dimKrylov << std::endl;
      this->SetSolver(defaultSolver);
    }
  }
  
  const auto sz = fStructMatrix->EquationFilter().NActiveEquations();
  auto &eigSolver = this->EigenSolver<TVar>();
  if (eigSolver.IsGeneralised()) {
    
    
    auto &materialVec = fCompMesh->MaterialVec();
    for (auto &&item : materialVec) {
      auto mat = (item.second);
      auto eigmat = dynamic_cast<TPZMatGeneralisedEigenVal *>(mat);
      if (eigmat)
        eigmat->SetMatrixA();
      else {
        PZError << __PRETTY_FUNCTION__;
        PZError << "\nERROR: cannot solve generalised EVP if\n"
                << "the materials do not have the interface:\n"
                << "TPZMatGeneralisedEigenVal\n"
                << "Aborting...\n";
        DebugStop();
      }
    }
  }
  auto matA = eigSolver.MatrixA();
  if (matA && matA->Rows() == sz) {
    matA->Zero();
    fStructMatrix->Assemble(*matA, dummyRhs, fGuiInterface);
  } else {
    TPZAutoPointer<TPZMatrix<TVar>> mat =
      dynamic_cast<TPZMatrix<TVar>*>(
          fStructMatrix->CreateAssemble(dummyRhs, fGuiInterface));
    eigSolver.SetMatrixA(mat);
  }
  PardisoConfig(eigSolver.MatrixA());
  // fSolver->UpdateFrom(fSolver->MatrixA());
  if (eigSolver.IsGeneralised()) {
    auto &materialVec = fCompMesh->MaterialVec();
    for (auto &&item : materialVec) {
      auto mat = (item.second);
      auto eigmat = dynamic_cast<TPZMatGeneralisedEigenVal *>(mat);
      if (eigmat)
        eigmat->SetMatrixB();
      else {
        PZError << __PRETTY_FUNCTION__;
        PZError << "\nERROR: cannot solve generalised EVP if\n"
                << "the materials do not have the interface:\n"
                << "TPZMatGeneralisedEigenVal\n"
                << "Aborting...\n";
        DebugStop();
      }
    }
    auto matB = eigSolver.MatrixB();
    if (matB && matB->Rows() == sz) {
      matB->Zero();
      fStructMatrix->Assemble(*matB, dummyRhs, fGuiInterface);
    } else {
      TPZAutoPointer<TPZMatrix<TVar>> mat =
        dynamic_cast<TPZMatrix<TVar>*>(
          fStructMatrix->CreateAssemble(dummyRhs, fGuiInterface));
      eigSolver.SetMatrixB(mat);
    }
    PardisoConfig(eigSolver.MatrixB());
  }
}

void TPZEigenAnalysis::Solve()
{
  if(fSolType == EReal)
    SolveT<STATE>();
  else
    SolveT<CSTATE>();
}

template<class TVar>
void TPZEigenAnalysis::SolveT()
{
    const auto nEq = fCompMesh->NEquations();
    const auto nReducedEq = fStructMatrix->NReducedEquations();
    const bool isReduced = nReducedEq == nEq ? false : true;
    auto &eigSolver = this->EigenSolver<TVar>();
    //if there is no equation filter nReducedEq == numeq
    fEigenvalues.Resize(nReducedEq);
    if (ComputeEigenvectors()) {
      TPZFMatrix<CTVar> *eigvectors = &fEigenvectors;
      if(isReduced){
        eigvectors = new TPZFMatrix<CTVar>;
      }
      const auto nev = eigSolver.NEigenpairs();
      eigvectors->Redim(nReducedEq, nev);
      eigSolver.Solve(fEigenvalues, *eigvectors);
      if(isReduced){
        fEigenvectors.Redim(nEq,nev);
        fStructMatrix->EquationFilter().Scatter(*eigvectors,fEigenvectors);
        delete eigvectors;
      }
    } else {
      eigSolver.Solve(fEigenvalues);
    }
}

TPZFMatrix<CSTATE>
TPZEigenAnalysis::GetEigenvectors() const {
  if(!ComputeEigenvectors()){
    std::cout<<__PRETTY_FUNCTION__
             <<"\nWARNING:No eigenvectors were computed\n";
  }
  return fEigenvectors;
}

TPZVec<CSTATE>
TPZEigenAnalysis::GetEigenvalues() const
{
    return fEigenvalues;
}

int TPZEigenAnalysis::ClassId() const
{
  return Hash("TPZEigenAnalysis") ^
    TPZAnalysis::ClassId() << 1;
}
  
void TPZEigenAnalysis::Write(TPZStream &buf, int withclassid) const
{
  TPZAnalysis::Write(buf,withclassid);
  buf.Write(fEigenvalues);
  fEigenvectors.Write(buf,withclassid);
  buf.Write(fCalcVectors);
}

void TPZEigenAnalysis::Read(TPZStream &buf, void *context)
{
  TPZAnalysis::Read(buf,context);
  buf.Read(fEigenvalues);
  fEigenvectors.Read(buf,context);
  buf.Read(fCalcVectors);
}


#define INSTANTIATE_TEMPLATES(TVar)                                            \
  template TPZEigenSolver<TVar> &TPZEigenAnalysis::EigenSolver<TVar>();        \

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

template class TPZRestoreClass<TPZEigenAnalysis>;
