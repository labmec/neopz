#include "TPZEigenAnalysisBase.h"
#include "TPZEigenSolver.h"
#include "pzcmesh.h"

TPZEigenAnalysisBase::TPZEigenAnalysisBase()
  : TPZRegisterClassId(&TPZEigenAnalysisBase::ClassId),TPZAnalysis()
{

}

TPZEigenAnalysisBase::TPZEigenAnalysisBase(TPZCompMesh *mesh,
                                   const RenumType &rt, std::ostream &out)
  : TPZAnalysis(mesh, rt,out)
{

}

TPZEigenAnalysisBase::TPZEigenAnalysisBase(TPZAutoPointer<TPZCompMesh> mesh,
                                           const RenumType &rt,
                                   std::ostream &out)
  : TPZAnalysis(mesh, rt, out)
{

}

int TPZEigenAnalysisBase::ClassId() const
{
  return Hash("TPZEigenAnalysisBase") ^
    TPZAnalysis::ClassId() << 1;
}
  
void TPZEigenAnalysisBase::Write(TPZStream &buf, int withclassid) const
{
  TPZAnalysis::Write(buf,withclassid);
  buf.Write(fEigenvalues);
  fEigenvectors.Write(buf,withclassid);
  buf.Write(fCalcVectors);
}

void TPZEigenAnalysisBase::Read(TPZStream &buf, void *context)
{
  TPZAnalysis::Read(buf,context);
  buf.Read(fEigenvalues);
  fEigenvectors.Read(buf,context);
  buf.Read(fCalcVectors);
}

void TPZEigenAnalysisBase::Solve()
{
  if(fSolType == EReal)
    SolveT<STATE>();
  else
    SolveT<CSTATE>();
}

template<class TVar>
void TPZEigenAnalysisBase::SolveT()
{
    const auto nEq = fCompMesh->NEquations();
    const auto nReducedEq = fStructMatrix->NReducedEquations();
    const bool isReduced = nReducedEq == nEq ? false : true;
    auto eigSolver = dynamic_cast<TPZEigenSolver<TVar>*>(this->Solver());
    if(!eigSolver){
        DebugStop();
    }
    const auto nev = eigSolver->NEigenpairs();
    fEigenvalues.Resize(nev);
    if (ComputeEigenvectors()) {
      TPZFMatrix<CTVar> *eigvectors = &fEigenvectors;
      if(isReduced){
        eigvectors = new TPZFMatrix<CTVar>;
      }
      //if there is no equation filter nReducedEq == numeq
      eigvectors->Redim(nReducedEq, nev);
      eigSolver->Solve(fEigenvalues, *eigvectors);
      if(isReduced){
        fEigenvectors.Redim(nEq,nev);
        fStructMatrix->EquationFilter().Scatter(*eigvectors,fEigenvectors);
        delete eigvectors;
      }
    } else {
      eigSolver->Solve(fEigenvalues);
    }
}