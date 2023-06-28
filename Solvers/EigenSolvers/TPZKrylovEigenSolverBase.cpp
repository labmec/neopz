#include "TPZKrylovEigenSolverBase.h"
#include "TPZLapackEigenSolver.h"
#include "TPZSimpleTimer.h"


template<class TVar>
bool TPZKrylovEigenSolverBase<TVar>::ArnoldiIteration(
  TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
  TPZFMatrix<TVar> &H, TVar &beta)
{

  if(KrylovDim() < 2){
    fKrylovDim = 10;
  }
  const int64_t nRows = this->SystemSize();
  const int n = std::min(fKrylovDim,nRows);


  if(Q.size() != n){
    Q.Resize(n,nullptr);
  }

  /*Q might already have a few computed vectors from previous iterations*/
  int first_k = 0;
  while(first_k < n && Q[first_k]){first_k++;}
  if(first_k == n){
    DebugStop();
  }

  if(H.Rows() != n || H.Cols() != n){
    H.Redim(n,n);
  }else{
    for(int i = first_k; i < n; i++){
      auto *hptr = &H.g(0,i);
      for(int row = 0; row < n; row++){
        *hptr++ = 0;
      }
    }
  }
  
  for(int i = first_k; i < n; i++) Q[i]= new TPZFMatrix<TVar>;
  
  /*see Chapter 2 of slepc manual(EPS) or search for Arnoldi Iteration*/

  //initializing first vector
  *(Q[first_k]) = fKrylovVector * (TVar)(1./Norm(fKrylovVector));
  
  TPZSimpleTimer arnoldiIteration("ArnoldiIteration");
  const auto tol = std::numeric_limits<RTVar>::epsilon();

  TPZFMatrix<TVar> w(nRows,1,0.);

  for(auto k = first_k+1; k < n+1; k++){
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
    }else{
      beta = normW;
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
  
  const int &krylovDim = this->KrylovDim();
  const auto size = this->SystemSize();
  TPZManVector<TPZAutoPointer<TPZFMatrix<TVar>>,20> qVecs(krylovDim, nullptr);
  TPZFNMatrix<400,TVar> h(krylovDim,krylovDim,0.);

  //this will be used as the first krylov vector
  if(fKrylovVector.Rows() != size || fKrylovVector.Cols() != 1){
    fKrylovVector.AutoFill(size,1,SymProp::NonSym);
  }

  auto myself = dynamic_cast<TPZEigenSolver<TVar>*>(this);
  if(!myself){
    DebugStop();
  }
  const int n = myself->NEigenpairs(); 
  bool all_converged = false;

  TVar beta{0};
  int lapackres{-1};
  int n_iter{0};
  const int max_iter = this->GetMaxIterations();
  const auto tol = this->Tolerance();
  TPZFNMatrix<400,CTVar> lapackEV(krylovDim,krylovDim,0.);
  TPZManVector<int,20> indices;

  std::cout<<"max it: "<<max_iter<<" tol: "<<tol<<std::endl;
  while(!all_converged && n_iter < max_iter){
    auto success = this->ArnoldiIteration(qVecs,h, beta);
    if(!success){
      return -1;
    }

    lapackEV.Zero();
  

    lapackres = [&h,&w,&lapackEV]()
    {
      TPZSimpleTimer lapacktimer("Hessenberg EVP");
      TPZLapackEigenSolver<TVar> lapack;
      return lapack.SolveHessenbergEigenProblem(h, w, lapackEV);
    }();
    if(lapackres) return lapackres;

    this->TransformEigenvalues(w);

    indices.Fill(0);
    myself->SortEigenvalues(w,indices);

    {
      std::set<int> converged_pairs;
      //minimum residual of non converged pairs
      STATE min_res{1e15};
      int min_index{0};
      for(int i = 0; i < n; i++){
        auto il = indices[i];
        const auto res = std::abs(beta*lapackEV.Get(krylovDim-1,il));
        if(res < tol){converged_pairs.insert(il);}
        else{
          if(res < min_res){
            min_res = res;
          }
        }
      }
      const int nconv = converged_pairs.size();
      std::cout<<"\riter "<<n_iter<<": "<<nconv<<" converged eigenpairs"
               <<" min res of non conv pair: "<<min_res<<std::flush;
      n_iter++;
      if(nconv >= n || n_iter >= max_iter){all_converged=true; break;}
      //now we compute the initial vector of the next arnoldi iteration
      //we choose it as the dominant eigenvector of H
      fKrylovVector.Zero();
      if constexpr (std::is_same_v<TVar,CTVar>){
        for (int j = 0; j < krylovDim; j++){//which vector from Q
          TVar *ev = &fKrylovVector.g(0,0);
          const auto lev = lapackEV(j,min_index);
          TVar *q = &qVecs[j]->g(0,0);
          /*
            The following loop computes
            eigenVectors(k,i) += lev * qvec.GetVal(k,0),
            with const auto qvec = *qVecs[j],
            but twice as fast
          */
          for(int k = 0; k < size; k++)
            *ev++ += lev * *q++;
        }
      }else{
        fKrylovVector.AutoFill(size,1,SymProp::NonSym);
      }
      //now we lock the converged eigenpairs
      TPZManVector<TPZAutoPointer<TPZFMatrix<TVar>>,20> qcp = qVecs;
      int count = 0;
      auto hcp = h;
      h.Zero();
      for(auto conv : converged_pairs){
        qVecs[count] = qcp[conv];
        int ccount = 0;
        for(auto cconv : converged_pairs){
          const auto val = hcp.Get(conv,cconv);
          h.Put(count,ccount,val);
          ccount++;
        }
        count++;
      }
      
      //now we erase all other vectors
      for(int i = nconv; i < krylovDim; i++){
        qVecs[i] = nullptr;
      }
    }
  }

  std::cout<<"\rfinished!"<<std::endl;
  for(int i = 0; i < n; i++){
    auto il = indices[i];
    const auto res = std::abs(beta*lapackEV.Get(krylovDim-1,il));
    std::cout<<"w :"<<w[i]<<" res "<<res<<std::endl;
  }
  
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