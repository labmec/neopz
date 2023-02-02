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
int TPZKrylovEigenSolver<TVar>::SolveImpl(TPZVec<CTVar> &w,
                                          TPZFMatrix<CTVar> &eigenVectors,
                                          bool computeVectors)
{
#ifdef USING_MKL
  auto PardisoSetup = [](auto pardiso_control, auto sys, auto prop, auto str){
    pardiso_control->SetStructure(str);
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
#else
    DebugStop();
#endif
  if(this->NEigenpairs() < 1) SetNEigenpairs(1);
#ifndef USING_LAPACK
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
#endif
  TPZSimpleTimer total("Arnoldi Solver",true);

  
  const int64_t nRows = this->MatrixA()->Rows();
  
  if(fUserTarget) AdjustTargetST();
  
  TPZAutoPointer<TPZMatrix<TVar>> arnoldiMat{nullptr};
  auto st = this->SpectralTransform();
  if(st){
    TPZSimpleTimer calcMat("ST Calculating matrix",true);
#ifdef USING_MKL
    auto prdscfg = this->GetPardisoControlA();
    if(prdscfg && !prdscfg->HasCustomSettings()){
      const auto str =
        this->MatrixA()->IsSymmetric() ?
        TPZPardisoSolver<TVar>::MStructure::ESymmetric:
        TPZPardisoSolver<TVar>::MStructure::ENonSymmetric;
      const auto sys =
        this->MatrixA()->IsSymmetric() ? 
        TPZPardisoSolver<TVar>::MSystemType::ESymmetric:
        TPZPardisoSolver<TVar>::MSystemType::ENonSymmetric;
      const auto prop =
        this->MatrixA()->IsDefPositive() ?
        TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
        TPZPardisoSolver<TVar>::MProperty::EIndefinite;
      PardisoSetup(prdscfg,sys,prop,str); 
    }
#else
      DebugStop();
#endif
    if(this->IsGeneralised())
      arnoldiMat = st->CalcMatrix(this->MatrixA(),this->MatrixB());
    else
      arnoldiMat = st->CalcMatrix(this->MatrixA());
  }else{
    arnoldiMat = this->MatrixA();
    if(this->IsGeneralised()){
      TPZSimpleTimer binvert("invert B mat",true);
#ifdef USING_MKL
      auto prdscfg = this->GetPardisoControlB();
      if(prdscfg && !prdscfg->HasCustomSettings()){
        const auto str =
        this->MatrixB()->IsSymmetric() ?
        TPZPardisoSolver<TVar>::MStructure::ESymmetric:
        TPZPardisoSolver<TVar>::MStructure::ENonSymmetric;
      const auto sys =
        this->MatrixB()->IsSymmetric() ? 
        TPZPardisoSolver<TVar>::MSystemType::ESymmetric:
        TPZPardisoSolver<TVar>::MSystemType::ENonSymmetric;
      const auto prop =
        this->MatrixB()->IsDefPositive() ?
        TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
        TPZPardisoSolver<TVar>::MProperty::EIndefinite;
        PardisoSetup(prdscfg,sys,prop,str); 
      }
#else
        DebugStop();
#endif
      if (this->MatrixB()->IsSymmetric()) this->MatrixB()->Decompose(ELDLt);
      else this->MatrixB()->Decompose(ELU);
    }
  }
  
  const int &n = this->NEigenpairs();
  if(KrylovDim() == -1){
    SetKrylovDim(10*n);
  }
  const int &krylovDim = KrylovDim();
  TPZManVector<TPZAutoPointer<TPZFMatrix<TVar>>,20> qVecs;
  TPZFNMatrix<400,TVar> h(krylovDim,krylovDim,0.);

  
  auto success = ArnoldiIteration(*arnoldiMat,qVecs,h);
  if(!success){
    return -1;
  }
  TPZFNMatrix<400,CTVar> lapackEV(n,n,0.);
  

  auto lapackres = [&h,&w,&lapackEV]()
  {
    TPZSimpleTimer lapacktimer("Hessenberg EVP",true);
    TPZLapackEigenSolver<TVar> lapack;
    return lapack.SolveHessenbergEigenProblem(h, w, lapackEV);
  }();
  if(lapackres) return lapackres;

  if(st) st->TransformEigenvalues(w);

  TPZManVector<int,20> indices;
  
  this->SortEigenvalues(w,indices);
  
  if(!computeVectors) return lapackres;
  eigenVectors.Redim(nRows,n);
  {
    TPZSimpleTimer evTimer("Computing eigenvectors",true);
    for (int i = 0; i< n; i++){//which eigenvector from A
      auto il = indices[i];
      for (int j = 0; j < krylovDim; j++){//which vector from Q
        CTVar *ev = &eigenVectors.g(0,i);
        const auto lev = lapackEV(j,il);
        TVar *q = &qVecs[j]->g(0,0);
        /*
          The following loop computes
          eigenVectors(k,i) += lev * qvec.GetVal(k,0),
          with const auto qvec = *qVecs[j],
          but twice as fast
        */
        for(int k = 0; k < nRows; k++)
          *ev++ += lev * *q++;
      }
    }
  }
  
  return lapackres;
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors)
{
  return SolveImpl(w,eigenVectors,true);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> eigenVectors;
  return SolveImpl(w,eigenVectors,false);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                 TPZFMatrix<CTVar> &eigenVectors)
{
  return SolveImpl(w,eigenVectors,true);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> eigenVectors;
  return SolveImpl(w, eigenVectors,false);
}


template<class TVar>
bool TPZKrylovEigenSolver<TVar>::ArnoldiIteration(
  const TPZMatrix<TVar> &A,
  TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
  TPZFMatrix<TVar> &H)
{

  if(KrylovDim() < 2){
    fKrylovDim = 10;
  }
  const int64_t nRows = A.Rows();
  const TPZMatrix<TVar> &B = this->fMatrixB.operator*();
  const int n = std::min(fKrylovDim,nRows);
  std::cout<<"Calculating Krylov subspace of dimension "<<n<<'\n';
  H.Redim(n,n);
  Q.Resize(n, nullptr);

  if(fKrylovVector.Rows() != nRows || fKrylovVector.Cols() != 1){
    fKrylovVector.AutoFill(nRows,1,0);
  }
  
  for(int i = 0; i < n; i++) Q[i]= new TPZFMatrix<TVar>;

  //deciding whether to multiply by b before, after and dont multiply at all
  enum class EWhichB{ENoB, EBBefore, EBAfter};

  EWhichB whichB = [this](){
    if(this->fIsGeneralised){
      auto st = this->SpectralTransform().operator->();
      auto stshiftinvert = dynamic_cast<TPZSTShiftAndInvert<TVar>*>(st);
      if(stshiftinvert) return EWhichB::EBBefore;
      else return EWhichB::EBAfter;
    }
    return EWhichB::ENoB;
  }();
  
  /*see Chapter 2 of slepc manual(EPS) or search for Arnoldi Iteration*/

  //initializing first vector
  *(Q[0]) = fKrylovVector * (TVar)(1./Norm(fKrylovVector));
  
  TPZSimpleTimer arnoldiIteration("ArnoldiIteration",true);
  const auto &tol = Tolerance();

  for(auto k = 1; k < n+1; k++){
    // TPZSimpleTimer arnoldiStep("step"+std::to_string(k),true);


    //let us generate a first guess for w: w = A.q_{k-1}
    TPZFMatrix<TVar> w = [&A,&B,&Q,k,whichB]()
    {
      // TPZSimpleTimer matMult("matmult",true);
      switch(whichB){
      case EWhichB::ENoB:
        return  A  * *(Q[k-1]);
      case EWhichB::EBBefore:
        return A * (B * *(Q[k-1]));
      case EWhichB::EBAfter:
        return B * (A * *(Q[k-1]));
      }
      unreachable();
    }();

    RTVar normW{1};
    bool success = false;
    /** after orthogonalising w.r.t. previous vectors (gram-schmidt)
        we will then have w_k = Av_k - sum_j^k (h_{jk} v_j)
    */
    while(!success){
      for(auto j = k-1; j >= 0; j--){
        const auto& qj = *(Q[j]);
        const auto dotqj = Dot(w,qj);
        H.PutVal(j,k-1,dotqj);
        w -= qj * dotqj;
      }

      normW = Norm(w);
      if(normW > tol || k == n){
        success = true;
      }
      else{
        //generate random unit vector and try again
        w.AutoFill(nRows,1,0);
        w *= 1/Norm(w);
      }
    }

    
    if(k < n){
      H.PutVal(k,k-1,normW);
      w *= (TVar)1./normW;
      (*(Q[k])) = std::move(w);
    }

    
  }//for k

  return true;
}

template class TPZKrylovEigenSolver<float>;
template class TPZKrylovEigenSolver<double>;
template class TPZKrylovEigenSolver<std::complex<float>>;
template class TPZKrylovEigenSolver<std::complex<double>>;
