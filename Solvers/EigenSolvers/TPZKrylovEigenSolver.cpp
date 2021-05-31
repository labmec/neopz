#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "TPZSimpleTimer.h"


#include <numeric>

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveImpl(TPZVec<CTVar> &w,
                                          TPZFMatrix<CTVar> &eigenVectors,
                                          bool computeVectors)
{
  
#ifndef USING_LAPACK
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
#endif
  TPZSimpleTimer total("ArnoldiSolver");
  auto &matA = this->fMatrixA.operator*();
  auto &matB = this->fMatrixB.operator*();
  const int nRows = matA.Rows();

  TPZAutoPointer<TPZMatrix<TVar>> arnoldiMat{nullptr};
  auto st = this->SpectralTransform();
  if(st){
    TPZSimpleTimer calcMat("ST calc mat");
    if(this->IsGeneralised())
      arnoldiMat = st->CalcMatrix(matA,matB);
    else
      arnoldiMat = st->CalcMatrix(matA);
  }else{
    if(this->IsGeneralised()){
      arnoldiMat = matA.Clone();
      TPZSimpleTimer binvert("invert B mat");
      if (matB.IsSymmetric()) matB.Decompose_LDLt();
      else matB.Decompose_LU();
    }
    else
      arnoldiMat = matA.Clone();
  }
  
  const int &n = NEigenpairs();
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

  const auto eigOrder = EigenSorting();
  auto sortFunc = [eigOrder](const CTVar a, const CTVar b){
    switch(eigOrder){
    case TPZEigenSort::EAbsAscending: return fabs(a) < fabs(b);
    case TPZEigenSort::EAbsDescending: return fabs(a) > fabs(b);
    case TPZEigenSort::ERealAscending: return a.real() < b.real();
    case TPZEigenSort::ERealDescending: return a.real() > b.real();
    case TPZEigenSort::EImagAscending: return a.imag() < b.imag();
    case TPZEigenSort::EImagDescending: return a.imag() > b.imag();
    }
    unreachable();
  };
  //sorting eigenvalues
  TPZVec<int> indices(krylovDim);
  std::iota(indices.begin(), indices.end(), 0); // Initializing
  std::stable_sort(
      indices.begin(), indices.end(),
      [&w, &sortFunc](int i, int j) { return sortFunc(w[i], w[j]); });
  std::stable_sort(w.begin(), w.end(),
                   [&sortFunc](auto i, auto j) { return sortFunc(i, j); });

  w.Resize(n);

  // for (int i = 0; i < n; i++)
  //   w[i] = (TVar)1.0/w[i] + shift;
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
  for(auto k = 0; k < n; k++){
    // TPZSimpleTimer arnoldiStep("step"+std::to_string(k));
      
    TPZFMatrix<TVar> w = [&A,&B,&Q,k,whichB]()
    {
      // TPZSimpleTimer matMult("matmult");
      switch(whichB){
      case EWhichB::ENoB:
        return  A  * *(Q[k]);
      case EWhichB::EBBefore:
        return A * (B * *(Q[k]));
      case EWhichB::EBAfter:
        return B * (A * *(Q[k]));
      }
      unreachable();
    }();
    /** after orthogonalising w.r.t. previous vectors (gram-schmidt)
        we will then have w_k = Av_k - sum_j^k (h_{jk} v_j)
    */
    {
      //tests indicated better precision if the loop is done with decreasing j
      // TPZSimpleTimer orth("orthogonalising");
      for(auto j = k; j >=0; j--){
        const auto& qj = *(Q[j]);
        H.PutVal(j,k,Dot(w,qj));
        w -= qj * H.GetVal(j,k);
      }
    }
      
    const auto normW = Norm(w);
    if(k<n-1) {
      H(k+1,k) = normW;
      w *= (TVar)1./normW;
      (*(Q[k+1])) = std::move(w);
    }
    if (normW < tol){
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nERROR:could not create krylov subspace\n";
      return false;
    }
  }//for k
  
  return true;
}

template class TPZKrylovEigenSolver<float>;
template class TPZKrylovEigenSolver<double>;
template class TPZKrylovEigenSolver<std::complex<float>>;
template class TPZKrylovEigenSolver<std::complex<double>>;