#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "TPZPardisoSolver.h"
#include "TPZSimpleTimer.h"

#include <algorithm>

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetTarget(TVar target)
{
  TPZEigenSolver<TVar>::SetTarget(target);
  fUserTarget=true;
  AdjustTargetST();
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::AdjustTargetST()
{
  auto st =
    dynamic_cast<TPZSTShiftOrigin<TVar>*>(this->SpectralTransform().operator->());
  if(st){
    st->SetShift(this->Target());
  }
}



template<class TVar>
void TPZKrylovEigenSolver<TVar>::PreSolve()
{
#ifdef USING_MKL
  auto PardisoSetup = [](auto pardiso_control, auto sys, auto prop){
    pardiso_control->SetMatrixType(sys,prop);
    auto param = pardiso_control->GetParam();
      //fParam[0] No default values
    param[0] = 1;
    //param[1]  use Metis for the ordering
    param[1] = 2;
    /*param[3]  Preconditioned CGS/CG. 
      0 = // No iterative-direct algorithm
      10*L+K
      L = stoppping criterion: 10^-L
      K = 
      0: The factorization is always computed as required by phase
      1: CGS iteration replaces the computation of LU. 
      The preconditioner is LU that was computed at a previous step
      (the first step or last step with a failure) in a sequence of
      solutions needed for identical sparsity patterns.
      2: CGS iteration for symmetric positive definite matrices
      Replaces the computation of LLt. The preconditioner is LLT
      that was computed at a previous step
      (the first step or last step with a failure)
      in a sequence of solutions needed for identical sparsity patterns. 
    */
    //param[4]  No user fill-in reducing permutation
    param[3] = 0;
    param[4] = 0;

    //param[7]  Maximum number of iterative refinement steps that the solver performs when perturbed pivots are obtained during the numerical factorization. 
    param[7] = 8;
      	
    //param[8]  Tolerance level for the relative residual in the iterative refinement process. (10^-{param[8]})
    param[8] = 23;
    //param[9]  Perturb the pivot elements with 1E-param[9]
    param[9] = 23;
    //param[10]  Use nonsymmetric permutation and scaling MPS
    param[10] = 0;

      
    //param[12]  Maximum weighted matching algorithm is switched-off (default for symmetric).
    param[12] = 0;
    //param[26] Whether to check matrix data
    param[26] = 1;
    //param[59]  Do not use OOC
    param[59] = 0;
    pardiso_control->SetParam(param);
  };
#endif
  if(this->NEigenpairs() < 1) SetNEigenpairs(1);
  
  if(this->KrylovDim() < this->NEigenpairs()) {
    this->SetKrylovDim(this->NEigenpairs());
  }
#ifndef USING_LAPACK
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
#endif


  
  const int64_t nRows = this->MatrixA()->Rows();
  
  if(fUserTarget) AdjustTargetST();
  
  auto st = this->SpectralTransform();
  if(st){
    TPZSimpleTimer calcMat("ST Calculating matrix",true);
#ifdef USING_MKL
    auto prdscfg = this->GetPardisoControlA();
    if(prdscfg && !prdscfg->HasCustomSettings()){
      const auto sys = this->MatrixA()->GetSymmetry();
      const auto prop =
        this->MatrixA()->IsDefPositive() ?
        TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
        TPZPardisoSolver<TVar>::MProperty::EIndefinite;
      PardisoSetup(prdscfg,sys,prop); 
    }
#else
      DebugStop();
#endif
    if(this->IsGeneralised())
      st->CalcMatrix(this->MatrixA(),this->MatrixB());
    else
      st->CalcMatrix(this->MatrixA());
  }else{
    if(this->IsGeneralised()){
      TPZSimpleTimer binvert("invert B mat",true);
#ifdef USING_MKL
      auto prdscfg = this->GetPardisoControlB();
      if(prdscfg && !prdscfg->HasCustomSettings()){
      const auto sys = this->MatrixB()->GetSymmetry();
      const auto prop =
        this->MatrixB()->IsDefPositive() ?
        TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
        TPZPardisoSolver<TVar>::MProperty::EIndefinite;
        PardisoSetup(prdscfg,sys,prop); 
      }
#else
        DebugStop();
#endif
      const auto sp= this->MatrixB()->GetSymmetry();
      const bool use_lu = sp == SymProp::Herm || (!std::is_same_v<TVar,RTVar> && sp == SymProp::Sym);
      if (use_lu) this->MatrixB()->Decompose(ELU);
      else this->MatrixB()->Decompose(ELDLt);
    }
  }

}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::TransformEigenvalues(TPZVec<CTVar> &w)
{
  auto st = SpectralTransform();
  if(st) st->TransformEigenvalues(w);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors)
{
  return this->SolveImpl(w,eigenVectors,true);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> eigenVectors;
  return this->SolveImpl(w,eigenVectors,false);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                 TPZFMatrix<CTVar> &eigenVectors)
{
  return this->SolveImpl(w,eigenVectors,true);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> eigenVectors;
  return this->SolveImpl(w, eigenVectors,false);
}


//! Applies (maybe matrix-free) operator on a given vector
template<class TVar>
void TPZKrylovEigenSolver<TVar>::ApplyOperator(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &res) const
{
  auto shift_invert = TPZAutoPointerDynamicCast<TPZSTShiftAndInvert<TVar>>(this->fST);
  if(this->fIsGeneralised){
    if(shift_invert){
      this->fMatrixB->Multiply(x,res);
      auto dectype = this->fMatrixA->IsDecomposed();
      this->fMatrixA->SolveDirect(res,dectype);
    }else{
      this->fMatrixA->Multiply(x,res);
      auto dectype = this->fMatrixB->IsDecomposed();
      this->fMatrixB->SolveDirect(res,dectype);
    }
  }else{
    if(shift_invert){
      res = x;
      auto dectype = this->fMatrixA->IsDecomposed();
      this->fMatrixA->SolveDirect(res,dectype);
    }else{
      this->fMatrixA->Multiply(x,res);
    }
  }
}
//! System size (number of rows)
template<class TVar>
int64_t TPZKrylovEigenSolver<TVar>::NRows() const
{
  return this->fMatrixA->Rows();
}

template class TPZKrylovEigenSolver<float>;
template class TPZKrylovEigenSolver<double>;
template class TPZKrylovEigenSolver<std::complex<float>>;
template class TPZKrylovEigenSolver<std::complex<double>>;
