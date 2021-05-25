#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "TPZSimpleTimer.h"

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors)
{
#ifndef USING_LAPACK
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
#endif
  TPZSimpleTimer total("ArnoldiSolver");
  auto &matA = this->fMatrixA.operator*();
  const int nRows = matA.Rows();
  TPZAutoPointer<TPZMatrix<TVar>> invMat = matA.NewMatrix();
  TPZAutoPointer<TPZMatrix<TVar>> shiftedMat = matA.Clone();
  //calculating A-sigmaB
  const auto &shift = Shift();
  for(int i = 0; i < nRows; i++) shiftedMat->PutVal(i,i,shiftedMat->GetVal(i,i)-shift);
  // if (shiftedMat->IsSymmetric()) shiftedMat->Decompose_LDLt();
  // else shiftedMat->Decompose_LU();
  
  const int &n = NEigenpairs();
  const int &krylovDim = KrylovDim();
  TPZManVector<TPZAutoPointer<TPZFMatrix<TVar>>,20> qVecs;
  TPZFNMatrix<400,TVar> h(krylovDim,krylovDim,0.);

  
  auto success = ArnoldiIteration(*shiftedMat,qVecs,h);
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
  // for (int i = 0; i < n; i++)
  //   w[i] = (TVar)1.0/w[i] + shift;

    
  eigenVectors.Redim(nRows,n);
  {
    TPZSimpleTimer evTimer("Computing eigenvectors");
    for (int i = 0; i< n; i++){//which eigenvector from A
      for (int j = 0; j < krylovDim; j++){//which vector from Q
        CTVar *ev = &eigenVectors.g(0,i);
        const auto lev = lapackEV(j,i);
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
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w)
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: Not yet implemented.\nAborting...\n";
  DebugStop();
  return -1;
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                 TPZFMatrix<CTVar> &eigenVectors)
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: Not yet implemented.\nAborting...\n";
  DebugStop();
  return -1;
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w)
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: Not yet implemented.\nAborting...\n";
  DebugStop();
  return -1;
}

template<class TVar>
bool TPZKrylovEigenSolver<TVar>::ArnoldiIteration(
  const TPZMatrix<TVar> &A,
  TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
  TPZFMatrix<TVar> &H)
{
  const int nRows = A.Rows();
  const int n = std::min(fKrylovDim,nRows);
  std::cout<<"Calculating Krylov subspace of dimension "<<n<<'\n';
  H.Redim(n,n);
  Q.Resize(n, nullptr);

  if(fKrylovVector.Rows() != nRows || fKrylovVector.Cols() != 1){
    fKrylovVector.AutoFill(nRows,1,0);
  }
  
  for(int i = 0; i < n; i++) Q[i]= new TPZFMatrix<TVar>;


  /*see Chapter 2 of slepc manual(EPS) or search for Arnoldi Iteration*/

  //initializing first vector
  *(Q[0]) = fKrylovVector * (TVar)(1./Norm(fKrylovVector));
  
  TPZSimpleTimer arnoldiIteration("ArnoldiIteration");
  const auto &tol = Tolerance();
  for(auto k = 0; k < n; k++){
    // TPZSimpleTimer arnoldiStep("step"+std::to_string(k));
      
    TPZFMatrix<TVar> w = [&A,&Q,k]()
    {
      // TPZSimpleTimer matMult("matmult");
      return  A  * *(Q[k]);//at this point qk is a unit vector
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