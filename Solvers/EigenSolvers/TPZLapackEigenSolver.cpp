#include "TPZLapackEigenSolver.h"
#include "TPZLapack.h"
#include "Hash/TPZHash.h"
#include "pzsbndmat.h"

using namespace std::complex_literals;

template<class TVar>
TPZLapackEigenSolver<TVar>::TPZLapackEigenSolver()
{
#ifndef USING_LAPACK
  PZError<<"TPZLapackEigenSolver is only available if using LAPACK.\n";
  PZError<<"Aborting...\n";
  DebugStop();
#endif
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::ClassId() const
{
  return Hash("TPZLapackEigenSolver") ^
    TPZEigenSolver<TVar>::ClassId() << 1;
}

template<class TVar>
TPZLapackEigenSolver<TVar>* TPZLapackEigenSolver<TVar>::Clone() const
{
  return new TPZLapackEigenSolver<TVar>(*this);
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w,
                                                  TPZFMatrix<CTVar> &eigenVectors)
{
  TPZFMatrix<TVar> *Afull =
    dynamic_cast<TPZFMatrix<TVar>*>(this->fMatrixA.operator->());
  TPZSBMatrix<TVar> *Absym =
    dynamic_cast<TPZSBMatrix<TVar>*>(this->fMatrixA.operator->());
  if(Afull){
    return SolveEigenProblem(*Afull,w,eigenVectors);
  }
  else if (Absym){
    return SolveEigenProblem(*Absym,w,eigenVectors);
  }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
  }
  return 0;
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveEigenProblem(TPZVec <CTVar> &w)
{
  TPZFMatrix<TVar> *Afull =
    dynamic_cast<TPZFMatrix<TVar>*>(this->fMatrixA.operator->());
  TPZSBMatrix<TVar> *Absym =
    dynamic_cast<TPZSBMatrix<TVar>*>(this->fMatrixA.operator->());
  if(Afull){
    return SolveEigenProblem(*Afull,w);
  }
  else if (Absym){
    return SolveEigenProblem(*Absym,w);
  }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
  }
  return 0;
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveGeneralisedEigenProblem(
    TPZVec <CTVar> &w, TPZFMatrix <CTVar> &eigenVectors)
{
  TPZFMatrix<TVar> *Afull =
    dynamic_cast<TPZFMatrix<TVar>*>(this->fMatrixA.operator->());
  TPZFMatrix<TVar> *Bfull =
    dynamic_cast<TPZFMatrix<TVar>*>(this->fMatrixB.operator->());
  TPZSBMatrix<TVar> *Absym =
    dynamic_cast<TPZSBMatrix<TVar>*>(this->fMatrixA.operator->());
  TPZSBMatrix<TVar> *Bbsym =
    dynamic_cast<TPZSBMatrix<TVar>*>(this->fMatrixB.operator->());
  if(Afull && Bfull){
    return SolveGeneralisedEigenProblem(*Afull,*Bfull,w,eigenVectors);
  }
  else if (Absym && Bbsym){
    return SolveGeneralisedEigenProblem(*Absym,*Bbsym,w,eigenVectors);
  }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec <CTVar> &w){
  TPZFMatrix<TVar> *Afull =
    dynamic_cast<TPZFMatrix<TVar>*>(this->fMatrixA.operator->());
  TPZFMatrix<TVar> *Bfull =
    dynamic_cast<TPZFMatrix<TVar>*>(this->fMatrixB.operator->());
  TPZSBMatrix<TVar> *Absym =
    dynamic_cast<TPZSBMatrix<TVar>*>(this->fMatrixA.operator->());
  TPZSBMatrix<TVar> *Bbsym =
    dynamic_cast<TPZSBMatrix<TVar>*>(this->fMatrixB.operator->());
  if(Afull && Bfull){
    return SolveGeneralisedEigenProblem(*Afull,*Bfull,w);
  }
  else if (Absym && Bbsym){
    return SolveGeneralisedEigenProblem(*Absym,*Bbsym,w);
  }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
  }
  return 1;
}


/*******************
*    TPZFMATRIX    *
*******************/
template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveEigenProblem(TPZFMatrix<TVar> &A,
                                                  TPZVec <CTVar> &eigenValues,
                                                  TPZFMatrix <CTVar> &eigenVectors,
                                                  bool calcVectors){
  if (A.Rows() != A.Cols()) {
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported dimensions for matrix A\nAborting...\n";
    DebugStop();
  }
  char jobvl[] = "N";
  char jobvr[] = "V";
  if(!calcVectors) jobvr[0]='N';
  
  TPZFMatrix<TVar> VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
  int dim = A.Rows();
  TVar testwork;
  int lwork = 10+20*dim;
  int info;

  TPZVec<TVar> work(lwork);

  eigenValues.Resize(dim,0.);
  if(calcVectors) eigenVectors.Redim(dim,dim);
#ifdef USING_LAPACK
  if constexpr(std::is_same_v<RTVar,TVar>){//real types
    TPZVec<TVar> realeigen(dim,0.);
    TPZVec<TVar> imageigen(dim,0.);
    if constexpr (std::is_same_v<TVar,float>){
      sgeev_(jobvl, jobvr, &dim, A.fElem, &dim, &realeigen[0],
             &imageigen[0], VL.fElem, &dim,
             VR.fElem, &dim, &work[0], &lwork, &info);
    }else if constexpr (std::is_same_v<TVar,double>){
      dgeev_(jobvl, jobvr, &dim, A.fElem, &dim, &realeigen[0],
             &imageigen[0], VL.fElem, &dim,
             VR.fElem, &dim, &work[0], &lwork, &info);
    }else{
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nERROR:Unsupported type\nAborting...\n";
      DebugStop();
    }
    for(int i = 0 ; i < dim ; i ++){
      eigenValues[i] = realeigen[i] + (CTVar)1i*imageigen[i];
    }
    if(calcVectors){
      for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
          for( int iV = 0 ; iV < dim ; iV++ ){
            eigenVectors(iV,i) = VR(iV,i);
          }
        }
        else{
          for( int iV = 0 ; iV < dim ; iV++ ){
            eigenVectors(iV,i) = VR(iV,i) + (CTVar)1i * VR(iV,i+1) ;
            eigenVectors(iV,i + 1) = VR(iV,i) - (CTVar)1i * VR(iV,i+1) ;
          }
          i++;
        }
      }
    }
  }else{
    TPZVec< RTVar > rwork( 2 * dim);
    if constexpr (std::is_same_v<TVar,std::complex<float>>){
      cgeev_(jobvl, jobvr, &dim, (varfloatcomplex*)A.fElem, &dim,
             (varfloatcomplex*)&eigenValues[0],
             (varfloatcomplex*)VL.fElem, &dim,
             (varfloatcomplex*)eigenVectors.fElem, &dim,
             (varfloatcomplex*)&work[0], &lwork, &rwork[0], &info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zgeev_(jobvl, jobvr, &dim, (vardoublecomplex*)A.fElem, &dim,
             (vardoublecomplex*)&eigenValues[0],
             (vardoublecomplex*)VL.fElem, &dim,
             (vardoublecomplex*)eigenVectors.fElem, &dim,
             (vardoublecomplex*)&work[0], &lwork, &rwork[0], &info);
    }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
    }
  }
#endif
  if (info != 0) {
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:LAPACK call returned with info: "<<info<<"\n";
    PZError<<"Aborting...\n";
    DebugStop();
  }

  return info;
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveEigenProblem(TPZFMatrix<TVar> &A,
                                                  TPZVec <CTVar> &eigenValues){
  TPZFNMatrix<1,CTVar> dummy(1,1,0.);
  constexpr bool calcVectors{false};
  return SolveEigenProblem(A,eigenValues,dummy,calcVectors);
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveGeneralisedEigenProblem(
    TPZFMatrix<TVar> &A,
    TPZFMatrix< TVar> &B ,
    TPZVec <CTVar> &eigenValues,
    TPZFMatrix <CTVar> &eigenVectors,
    bool calcVectors)
{
  if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
  {
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"Incompatible dimensions\nAborting...\n";
    DebugStop();
  }

  char jobvl[] = "N", jobvr[] = "V";
  if(!calcVectors) jobvr[0] = 'N';
  TPZFMatrix< TVar> VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
  int dim = A.Rows();
  TVar testwork;
  int lwork = 10+20*dim;
  int info;

  TPZVec<TVar> beta(dim);
  TPZVec<TVar> work(lwork);

#ifdef USING_LAPACK
  
  if(calcVectors) eigenVectors.Redim(dim,dim);
  eigenValues.Resize(dim,0.);
  
  if constexpr(std::is_same_v<TVar,RTVar>){//real types
    
    TPZVec<TVar> realeigen(dim,0.);
    TPZVec<TVar> imageigen(dim,0.);
    if constexpr (std::is_same_v<TVar,float>){
      sggev_(jobvl, jobvr, &dim, A.fElem, &dim , B.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);
    }else if constexpr (std::is_same_v<TVar,double>){
      dggev_(jobvl, jobvr, &dim, A.fElem, &dim , B.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);
    }

    for(int i = 0 ; i < dim ; i ++){
      if( IsZero(beta[i])){
        DebugStop(); //fran: i really dont know what to do with this result
      }
      else{
        eigenValues[i] = (realeigen[i] + (CTVar)1i*imageigen[i]) / beta[i];
      }
    }
    if(calcVectors){
      for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
          for( int iV = 0 ; iV < dim ; iV++ ){
            eigenVectors(iV,i) = VR(iV,i);
          }
        }
        else{
          for( int iV = 0 ; iV < dim ; iV++ ){
            eigenVectors(iV,i) = VR(iV,i) + (CTVar)1i * VR(iV,i+1) ;
            eigenVectors(iV,i + 1) = VR(iV,i) - (CTVar)1i * VR(iV,i+1) ;
          }
          i++;
        }
      }
    }
  }else if constexpr (std::is_same_v<TVar,CTVar>){//complex types
    TPZVec<TVar> eigen(dim,0.);
    TPZVec<RTVar> rwork( 8 * dim );
    if constexpr (std::is_same_v<TVar,std::complex<float>>){
      cggev_(jobvl, jobvr, &dim,
             (varfloatcomplex*)A.fElem, &dim,
             (varfloatcomplex*)B.fElem, &dim,
             (varfloatcomplex*)&eigen[0], (varfloatcomplex*)&beta[0],
             (varfloatcomplex*)VL.fElem, &dim,
             (varfloatcomplex*)eigenVectors.fElem, &dim,
             (varfloatcomplex*)&work[0], &lwork, &rwork[0],&info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zggev_(jobvl, jobvr, &dim,
             (vardoublecomplex*)A.fElem, &dim,
             (vardoublecomplex*)B.fElem, &dim,
             (vardoublecomplex*)&eigen[0], (vardoublecomplex*)&beta[0],
             (vardoublecomplex*)VL.fElem, &dim,
             (vardoublecomplex*)eigenVectors.fElem, &dim,
             (vardoublecomplex*)&work[0], &lwork, &rwork[0],&info);
    }

    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenValues[i] = eigen[i] / beta[i];
        }
    }
  }else{
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nERROR:Unsupported type\nAborting...\n";
      DebugStop();
    }
#endif

  if (info != 0) {
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:LAPACK call returned with info: "<<info<<"\n";
    PZError<<"Aborting...\n";
    DebugStop();
  }
  return info;
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveGeneralisedEigenProblem(
    TPZFMatrix<TVar> &A,
    TPZFMatrix<TVar> &B ,
    TPZVec<CTVar> &eigenValues)
{
  TPZFNMatrix<1,CTVar> dummy(1,1,0.);
  constexpr bool calcVectors{false};
  return SolveGeneralisedEigenProblem(A,B,eigenValues,dummy,calcVectors);
}

/*******************
*    TPZSBMATRIX    *
*******************/
template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveEigenProblem(TPZSBMatrix<TVar> &A,
                                                  TPZVec <CTVar> &eigenValues,
                                                  TPZFMatrix <CTVar> &eigenVectors,
                                                  bool calcVectors){
  if (A.Rows() != A.Cols())
  {
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"Incompatible dimensions\nAborting...\n";
    DebugStop();
  }

  char jobz = calcVectors ? 'V' : 'N'; //compute eigenvectors
  char uplo = 'U';//assume upper triangular
  int n = A.Dim();
  int kd = A.fBand;
  int ldab = A.fBand + 1;
  int ldbb = A.fBand + 1;
  TPZVec<RTVar> w(0,0.);
  w.Resize(n);
  int ldz = n;

  int info = -666;
#ifdef USING_LAPACK
  eigenValues.Resize(n);
  if(calcVectors) eigenVectors.Redim(n, n);
  
  if constexpr(std::is_same_v<TVar,RTVar>){
    TPZVec<TVar> work(3*n);
    TPZFMatrix<TVar> z(n,n);
    if constexpr (std::is_same_v<TVar,float>){
    
      ssbev_(&jobz, &uplo, &n, &kd, A.fDiag.begin(), &ldab, w.begin(),
             &z(0,0), &ldz, work.begin(), &info);
    }else if constexpr (std::is_same_v<TVar,double>){
      TPZVec<TVar> work(3*n);
      dsbev_(&jobz, &uplo, &n, &kd, A.fDiag.begin(), &ldab, w.begin(),
             &z(0,0), &ldz, work.begin(), &info);
    }
    if(calcVectors){
      for (int iVec = 0 ; iVec < n; iVec++) {
        for (int iCol = 0; iCol < n; iCol++) {
          eigenVectors( iVec , iCol) = z(iVec,iCol);
        }
      }
    }
  }else if constexpr(std::is_same_v<TVar,CTVar>){
    TPZVec<TVar> work(n);
    TPZVec<RTVar> rwork(3*n);
    if constexpr (std::is_same_v<TVar,std::complex<float>>){
      chbev_(&jobz, &uplo, &n, &kd,
             (varfloatcomplex*)A.fDiag.begin(), &ldab, w.begin(),
             (varfloatcomplex*)&eigenVectors(0,0), &ldz,
             (varfloatcomplex*)work.begin(), rwork.begin(), &info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zhbev_(&jobz, &uplo, &n, &kd,
             (vardoublecomplex*)A.fDiag.begin(), &ldab, w.begin(),
             (vardoublecomplex*)&eigenVectors(0,0), &ldz,
             (vardoublecomplex*)work.begin(), rwork.begin(), &info);
    }
  }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
  }
#endif

  if(info != 0){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:LAPACK call returned with info: "<<info<<"\n";
    PZError<<"Aborting...\n";
    DebugStop();
  }

  for(int i = 0 ; i < n ; i++){
    eigenValues[i] = w[i];
  }
  return info;
}
  
template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveEigenProblem(TPZSBMatrix<TVar> &A,
                                                  TPZVec <CTVar> &eigenValues){
  TPZFNMatrix<1,CTVar> dummy(1,1,0.);
  constexpr bool calcVectors{false};
  return SolveEigenProblem(A,eigenValues,dummy,calcVectors);
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveGeneralisedEigenProblem(
    TPZSBMatrix<TVar> &A, TPZSBMatrix< TVar> &B ,
    TPZVec <CTVar> &eigenValues, TPZFMatrix <CTVar> &eigenVectors,
    bool calcVectors)
{  
  if (  A.Rows() != B.Rows() && A.Cols() != B.Cols() )
  {
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported dimensions for matrix A\nAborting...\n";
    DebugStop();
  }


  char jobz = calcVectors? 'V' : 'N'; //Compute eigenvectors
  char uplo = 'U';//assume upper triangular
  int n = A.Dim();
  int ka = A.fBand;
  int kb = B.fBand;
  int ldab = A.fBand + 1;
  int ldbb = A.fBand + 1;
  TPZVec<RTVar> w(0,0.);
  w.Resize( n );
  int ldz = n;
  int info = -666;
#ifdef USING_LAPACK
  eigenValues.Resize(n);
  if(calcVectors) eigenVectors.Redim(n, n);
  
  if constexpr(std::is_same_v<TVar,RTVar>){
    TPZVec<TVar> work(3*n);
    TPZFMatrix<TVar> z(n,n);
    if constexpr (std::is_same_v<TVar,float>){
    
      ssbgv_(&jobz, &uplo, &n, &ka, &kb, A.fDiag.begin(),
             &ldab, B.fDiag.begin(), &ldbb, w.begin(),
             &z(0,0), &ldz, work.begin(), &info);
    }else if constexpr (std::is_same_v<TVar,double>){
      TPZVec<TVar> work(3*n);
      dsbgv_(&jobz, &uplo, &n, &ka, &kb, A.fDiag.begin(),
             &ldab, B.fDiag.begin(), &ldbb, w.begin(),
             &z(0,0), &ldz, work.begin(), &info);
    }
    if(calcVectors){
      for (int iVec = 0 ; iVec < n; iVec++) {
        for (int iCol = 0; iCol < n; iCol++) {
          eigenVectors( iVec , iCol) = z(iVec,iCol);
        }
      }
    }
  }else if constexpr(std::is_same_v<TVar,CTVar>){
    TPZVec<TVar> work(n);
    TPZVec<RTVar> rwork(3*n);
    if constexpr (std::is_same_v<TVar,std::complex<float>>){
      chbgv_(&jobz, &uplo, &n, &ka, &kb,
             (varfloatcomplex *)A.fDiag.begin(), &ldab,
             (varfloatcomplex *)B.fDiag.begin(), &ldbb, w.begin(),
             (varfloatcomplex *)&eigenVectors(0,0), &ldz,
             (varfloatcomplex *)work.begin(),rwork.begin(), &info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zhbgv_(&jobz, &uplo, &n, &ka, &kb,
             (vardoublecomplex *)A.fDiag.begin(), &ldab,
             (vardoublecomplex *)B.fDiag.begin(), &ldbb, w.begin(),
             (vardoublecomplex *)&eigenVectors(0,0), &ldz,
             (vardoublecomplex *)work.begin(),rwork.begin(), &info);
    }
  }else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:Unsupported type\nAborting...\n";
    DebugStop();
  }
#endif
  if(info != 0){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR:LAPACK call returned with info: "<<info<<"\n";
    PZError<<"Aborting...\n";
    DebugStop();
  }

  for(int i = 0 ; i < n ; i++){
    eigenValues[i] = w[i];
  }
  return info;
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveGeneralisedEigenProblem(
    TPZSBMatrix<TVar> &A,TPZSBMatrix< TVar> &B ,
    TPZVec <CTVar> &eigenValues)
{
  TPZFNMatrix<1,CTVar> dummy(1,1,0.);
  constexpr bool calcVectors{false};
  return SolveGeneralisedEigenProblem(A,B,eigenValues,dummy,calcVectors);
}


template class TPZLapackEigenSolver<std::complex<float>>;
template class TPZLapackEigenSolver<std::complex<double>>;

template class TPZLapackEigenSolver<float>;
template class TPZLapackEigenSolver<double>;