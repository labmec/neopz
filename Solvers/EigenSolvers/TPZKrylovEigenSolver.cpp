#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "TPZSimpleTimer.h"

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
  if(this->NEigenpairs() < 1) SetNEigenpairs(1);
#ifndef USING_LAPACK
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
#endif
  TPZSimpleTimer total("Arnoldi Solver");

  
  const int nRows = this->MatrixA()->Rows();
  
  if(fUserTarget) AdjustTargetST();
  
  TPZAutoPointer<TPZMatrix<TVar>> arnoldiMat{nullptr};
  auto st = this->SpectralTransform();
  if(st){
    TPZSimpleTimer calcMat("ST Calculating matrix");
    if(this->IsGeneralised())
      arnoldiMat = st->CalcMatrix(this->MatrixA(),this->MatrixB());
    else
      arnoldiMat = st->CalcMatrix(this->MatrixA());
  }else{
    arnoldiMat = this->MatrixA();
    if(this->IsGeneralised()){
      TPZSimpleTimer binvert("invert B mat");
      if (this->MatrixB()->IsSymmetric()) this->MatrixB()->Decompose_LDLt();
      else this->MatrixB()->Decompose_LU();
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
    TPZSimpleTimer lapacktimer("Hessenberg EVP");
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
    TPZSimpleTimer evTimer("Computing eigenvectors");
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
  const int nRows = A.Rows();
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
  
  TPZSimpleTimer arnoldiIteration("ArnoldiIteration");
  const auto &tol = Tolerance();

  for(auto k = 1; k < n+1; k++){
    // TPZSimpleTimer arnoldiStep("step"+std::to_string(k));


    //let us generate a first guess for w: w = A.q_{k-1}
    TPZFMatrix<TVar> w = [&A,&B,&Q,k,whichB]()
    {
      // TPZSimpleTimer matMult("matmult");
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