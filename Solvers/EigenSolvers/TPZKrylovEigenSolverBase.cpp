#include "TPZKrylovEigenSolverBase.h"
#include "TPZLapackEigenSolver.h"
#include "TPZSimpleTimer.h"


template<class TVar>
bool TPZKrylovEigenSolverBase<TVar>::ArnoldiIteration(
  TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
  TPZFMatrix<TVar> &H)
{

  if(KrylovDim() < 2){
    fKrylovDim = 10;
  }
  const int64_t nRows = this->SystemSize();
  const int n = std::min(fKrylovDim,nRows);
  std::cout<<"Calculating Krylov subspace of dimension "<<n<<'\n';


  if(Q.size() != n){
    Q.Resize(n,nullptr);
  }

  /*Q might already have a few computed vectors from previous iterations*/
  int first_k = 0;
  while(first_k < n && Q[first_k]){first_k++;}
  if(first_k == n){
    DebugStop();
  }
  H.Redim(n,n);
  Q.Resize(n, nullptr);

  if(fKrylovVector.Rows() != nRows || fKrylovVector.Cols() != 1){
    fKrylovVector.AutoFill(nRows,1,SymProp::NonSym);
  }
  
  for(int i = 0; i < n; i++) Q[i]= new TPZFMatrix<TVar>;
  
  /*see Chapter 2 of slepc manual(EPS) or search for Arnoldi Iteration*/

  //initializing first vector
  *(Q[0]) = fKrylovVector * (TVar)(1./Norm(fKrylovVector));
  
  TPZSimpleTimer arnoldiIteration("ArnoldiIteration");
  const auto &tol = Tolerance();

  TPZFMatrix<TVar> w(nRows,1,0.);

  for(auto k = 1; k < n+1; k++){
    // TPZSimpleTimer arnoldiStep("step"+std::to_string(k),true);

    //let us generate a first guess for w: w = A.q_{k-1}
    this->ApplyOperator(*Q[k-1],w);

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
        w.AutoFill(nRows,1,SymProp::NonSym);
        w *= 1/Norm(w);
      }
    }

    
    if(k < n){
      H.PutVal(k,k-1,normW);
      w *= (TVar)1./normW;
      (*(Q[k])) = w;
    }

    
  }//for k

  return true;
}

template<class TVar>
int TPZKrylovEigenSolverBase<TVar>::SolveImpl(TPZVec<CTVar> &w,
                                              TPZFMatrix<CTVar> &eigenVectors,
                                              bool computeVectors)
{
#ifndef USING_LAPACK
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
#endif
  TPZSimpleTimer total("KrylovEigenSolver::SolveImpl",true);

  this->PreSolve();
 
  if(this->KrylovDim() == -1){
    DebugStop();
  }
  
  int krylovDim = this->KrylovDim();
  const auto size = this->SystemSize();
  if(krylovDim > size){
    krylovDim=size;
    this->SetKrylovDim(size);
  }
  TPZManVector<TPZAutoPointer<TPZFMatrix<TVar>>,20> qVecs(krylovDim, nullptr);
  TPZFNMatrix<400,TVar> h(krylovDim,krylovDim,0.);

  if(fKrylovVector.Rows() != size || fKrylovVector.Cols() != 1){
    fKrylovVector.AutoFill(size,1,SymProp::NonSym);
  }
  //we initialise the first vector
  qVecs[0] = new TPZFMatrix<TVar>(fKrylovVector * (TVar)(1./Norm(fKrylovVector)));
  
  
  auto success = this->ArnoldiIteration(qVecs,h);
  if(!success){
    return -1;
  }


  auto myself = dynamic_cast<TPZEigenSolver<TVar>*>(this);
  if(!myself){
    DebugStop();
  }
  
  const int n = myself->NEigenpairs();
  
  TPZFNMatrix<400,CTVar> lapackEV(n,n,0.);
  

  auto lapackres = [&h,&w,&lapackEV]()
  {
    TPZSimpleTimer lapacktimer("Hessenberg EVP");
    TPZLapackEigenSolver<TVar> lapack;
    return lapack.SolveHessenbergEigenProblem(h, w, lapackEV);
  }();
  if(lapackres) return lapackres;

  this->TransformEigenvalues(w);

  TPZManVector<int,20> indices;


  
  myself->SortEigenvalues(w,indices);
  
  if(!computeVectors) return lapackres;
  const auto nRows = this->NRows();
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

template class TPZKrylovEigenSolverBase<float>;
template class TPZKrylovEigenSolverBase<double>;
template class TPZKrylovEigenSolverBase<std::complex<float>>;
template class TPZKrylovEigenSolverBase<std::complex<double>>;