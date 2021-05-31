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

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveHessenbergEigenProblem(TPZFMatrix<TVar> &A,
                                                            TPZVec<CTVar> &w,
                                                            TPZFMatrix<CTVar> &vecs)
{
  return SolveHessenbergEigenProblem(A,w,vecs,true);
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveHessenbergEigenProblem(TPZFMatrix<TVar> &A,
                                                            TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> vecs;
  return SolveHessenbergEigenProblem(A,w,vecs,false);
}

template<class TVar>
int TPZLapackEigenSolver<TVar>::SolveHessenbergEigenProblem(TPZFMatrix<TVar> &A,
                                                            TPZVec<CTVar> &w,
                                                            TPZFMatrix<CTVar> &vecs,
                                                            bool calcVecs)
{
  const int nrows = A.Rows();
  if(nrows != A.Cols()){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: mat is not a square matrix\n";
    return -1;
  }
#ifdef PZDEBUG
  for(auto i = 0; i < nrows; i++){
    for(int j = 0; j < i -1; j++)
      if(fabs(A.GetVal(i,j))!=0.){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: mat is not hessenberg\n";
        return -1;
      }
  }
#endif

  const char job = 'E';// compute eigenvalues only
  const char compz = 'N';//do not compute Schur vectors
  const int &n = nrows;//order of the matrix
  //ilo and ihi should be 1 and N if the matrix wasnt computed by ZGEBAL
  const int ilo = 1;
  const int ihi = n;
  auto H = A;//the matrix (is unspecified on exit)
  const int ldh = n;//leading dimension of the array
  //WR is set only for real types
  //WI is set only for real types
  //W is set only for complex types
  
  TPZVec<TVar> Z(nrows);//not referenced since compz=N
  const int ldz{1};//leading dimension of Z
  const int lwork{11*nrows};//dimension of work
  TPZVec<TVar> work(lwork);
  int info;

  w.Resize(n);

  TPZManVector<int,20> select(n,1);//which eigenvectors to compute in the next step
  TPZVec<TVar> wrVec, wiVec;
#ifdef USING_LAPACK
  if constexpr (std::is_same_v<TVar,RTVar>){
    wrVec.Resize(n);
    wiVec.Resize(n);
    if constexpr (std::is_same_v<TVar,float>){
      shseqr_(&job,&compz,&n,&ilo,&ihi,H.fElem,&ldh,&wrVec[0],&wiVec[0],&Z[0],&ldz,&work[0],&lwork,&info);
    }else if constexpr (std::is_same_v<TVar,double>){
      dhseqr_(&job,&compz,&n,&ilo,&ihi,H.fElem,&ldh,&wrVec[0],&wiVec[0],&Z[0],&ldz,&work[0],&lwork,&info);
    }
    for(int i = 0 ; i < n ; i ++){
      w[i] = wrVec[i] + (CTVar)1i*wiVec[i];
      if(!IsZero(wiVec[i])){
        select[i+1] = 0;//the other one is the complex conjugate
        w[i+1] = std::conj(w[i]);
        i++;
      }
    }
  }else{
    if constexpr (std::is_same_v<TVar,std::complex<float>>){
      chseqr_(&job,&compz,&n,&ilo,&ihi,
              (varfloatcomplex*)H.fElem,&ldh,
              (varfloatcomplex*)&w[0],
              (varfloatcomplex*)&Z[0],&ldz,
              (varfloatcomplex*)&work[0],&lwork,&info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zhseqr_(&job,&compz,&n,&ilo,&ihi,
              (vardoublecomplex*)H.fElem,&ldh,
              (vardoublecomplex*)&w[0],
              (vardoublecomplex*)&Z[0],&ldz,
              (vardoublecomplex*)&work[0],&lwork,&info);
    }
  }
#endif

  if(info != 0){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"Lapack returned with info : "<<info<<std::endl;
    return info;
  }

  if(!calcVecs) return info;
  
  const int mm = n;
  const char side{'R'};//compute right eigenvectors only
  const char eigsrc{'Q'};//eigenvalues were found from zhseqr
  const char initv{'N'};//no initial vectors
  //select has been set already
  //n has been set already
  H=A;
  //ldh has been set already
  //WR is set only for real types
  //WI is set only for real types
  //W is set only for complex types
  TPZVec<TVar> vl(mm);//not referenced
  const int ldvl{1};//will be ignored
  //vr will be the eigenvectors
  const int ldvr{n};//leading dimension of vr
  int m;
  //mm has already been set
  const int worksize = std::is_same_v<TVar,RTVar> ? n*n+2*n : n*n;
  work.Resize(worksize);
  vecs.Redim(n,n);
  TPZVec<int> ifaill(mm,-1);//will not be referenced
  TPZVec<int> ifailr(mm,-1);


#ifdef USING_LAPACK
  if constexpr(std::is_same_v<TVar,RTVar>){//real types
    TPZFMatrix<TVar> VR(ldvr,mm,-1);
    
    if constexpr(std::is_same_v<TVar,float>){
      shsein_(&side, &eigsrc, &initv, &select[0], &n,
              H.fElem, &ldh, &wrVec[0],&wiVec[0],
              &vl[0], &ldvl,VR.fElem, &ldvr, &mm, &m,
              &work[0],&ifaill[0], &ifailr[0], &info);
    }else if constexpr(std::is_same_v<TVar,double>){
      dhsein_(&side, &eigsrc, &initv, &select[0], &n,
              H.fElem, &ldh, &wrVec[0], &wiVec[0],
              &vl[0], &ldvl,VR.fElem, &ldvr, &mm, &m,
              &work[0],&ifaill[0], &ifailr[0], &info);
    }
    // std::cout<<"m "<<m<<" mm "<<mm<<std::endl;
    for(int i = 0 ; i < n ; i ++){
      for( int iV = 0 ; iV < n ; iV++ ){
        vecs(iV,i) = VR(iV,i);
      }
      if(i<n-1 && select[i+1] == 0){
        for( int iV = 0 ; iV < n ; iV++ ){
          vecs(iV,i + 1) = std::conj(vecs(iV,i));
        }
        i++;
      }
    }
  }else{//complex types
    TPZManVector<RTVar,20> rworkvec(n,0.);
    RTVar *rwork = &rworkvec[0];
    auto wCopy = w;
    if constexpr(std::is_same_v<TVar,std::complex<float>>){
      chsein_(&side,&eigsrc,&initv,&select[0],&n,
              (varfloatcomplex*)H.fElem,&ldh,
              (varfloatcomplex*)&wCopy[0],
              (varfloatcomplex*)&vl[0],&ldvl,
              (varfloatcomplex*)&vecs[0],&ldvr,&mm,&m,
              (varfloatcomplex*)&work[0],rwork,
              &ifaill[0],&ifailr[0],&info);
    }else if constexpr(std::is_same_v<TVar,std::complex<double>>){
      zhsein_(&side,&eigsrc,&initv,&select[0],&n,
              (vardoublecomplex*)H.fElem,&ldh,
              (vardoublecomplex*)&wCopy[0],
              (vardoublecomplex*)&vl[0],&ldvl,
              (vardoublecomplex*)&vecs[0],&ldvr,&mm,&m,
              (vardoublecomplex*)&work[0],rwork,
              &ifaill[0],&ifailr[0],&info);
    }
  }
#endif
  if(info != 0){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"Lapack returned with info : "<<info<<std::endl;
  }

  return info;
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
  TPZFMatrix <CTVar> eigenVectorsLapack;
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
  if(calcVectors) eigenVectorsLapack.Redim(dim,dim);
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
            eigenVectorsLapack(iV,i) = VR(iV,i);
          }
        }
        else{
          for( int iV = 0 ; iV < dim ; iV++ ){
            eigenVectorsLapack(iV,i) = VR(iV,i) + (CTVar)1i * VR(iV,i+1) ;
            eigenVectorsLapack(iV,i + 1) = VR(iV,i) - (CTVar)1i * VR(iV,i+1) ;
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
             (varfloatcomplex*)eigenVectorsLapack.fElem, &dim,
             (varfloatcomplex*)&work[0], &lwork, &rwork[0], &info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zgeev_(jobvl, jobvr, &dim, (vardoublecomplex*)A.fElem, &dim,
             (vardoublecomplex*)&eigenValues[0],
             (vardoublecomplex*)VL.fElem, &dim,
             (vardoublecomplex*)eigenVectorsLapack.fElem, &dim,
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
  const auto nev = this->NEigenpairs();
  if(nev > 0){
    TPZManVector<int,20> indices;
    this->SortEigenvalues(eigenValues,indices);
    if(calcVectors){
      eigenVectors.Resize(dim, nev);
      for (auto i = 0; i < nev; i++) {
        auto li = indices[i];
        for (auto x = 0; x < dim; x++)
          eigenVectors(x, i) = eigenVectorsLapack(x, li);
      }
    }
  }else{
    if(calcVectors) eigenVectors = std::move(eigenVectorsLapack);
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
  TPZFMatrix <CTVar> eigenVectorsLapack;
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
  
  if(calcVectors) eigenVectorsLapack.Redim(dim,dim);
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
            eigenVectorsLapack(iV,i) = VR(iV,i);
          }
        }
        else{
          for( int iV = 0 ; iV < dim ; iV++ ){
            eigenVectorsLapack(iV,i) = VR(iV,i) + (CTVar)1i * VR(iV,i+1) ;
            eigenVectorsLapack(iV,i + 1) = VR(iV,i) - (CTVar)1i * VR(iV,i+1) ;
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
             (varfloatcomplex*)eigenVectorsLapack.fElem, &dim,
             (varfloatcomplex*)&work[0], &lwork, &rwork[0],&info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zggev_(jobvl, jobvr, &dim,
             (vardoublecomplex*)A.fElem, &dim,
             (vardoublecomplex*)B.fElem, &dim,
             (vardoublecomplex*)&eigen[0], (vardoublecomplex*)&beta[0],
             (vardoublecomplex*)VL.fElem, &dim,
             (vardoublecomplex*)eigenVectorsLapack.fElem, &dim,
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

  const auto nev = this->NEigenpairs();
  if(nev > 0){
    TPZManVector<int,20> indices;
    this->SortEigenvalues(eigenValues,indices);
    if(calcVectors){
      eigenVectors.Resize(dim, nev);
      for (auto i = 0; i < nev; i++) {
        auto li = indices[i];
        for (auto x = 0; x < dim; x++)
          eigenVectors(x, i) = eigenVectorsLapack(x, li);
      }
    }
  }else{
    if(calcVectors) eigenVectors = std::move(eigenVectorsLapack);
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
  TPZFMatrix <CTVar> eigenVectorsLapack;
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
  if(calcVectors) eigenVectorsLapack.Redim(n, n);
  
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
          eigenVectorsLapack( iVec , iCol) = z(iVec,iCol);
        }
      }
    }
  }else if constexpr(std::is_same_v<TVar,CTVar>){
    TPZVec<TVar> work(n);
    TPZVec<RTVar> rwork(3*n);
    if constexpr (std::is_same_v<TVar,std::complex<float>>){
      chbev_(&jobz, &uplo, &n, &kd,
             (varfloatcomplex*)A.fDiag.begin(), &ldab, w.begin(),
             (varfloatcomplex*)&eigenVectorsLapack(0,0), &ldz,
             (varfloatcomplex*)work.begin(), rwork.begin(), &info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zhbev_(&jobz, &uplo, &n, &kd,
             (vardoublecomplex*)A.fDiag.begin(), &ldab, w.begin(),
             (vardoublecomplex*)&eigenVectorsLapack(0,0), &ldz,
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
  const auto nev = this->NEigenpairs();
  if(nev > 0){
    TPZManVector<int,20> indices;
    this->SortEigenvalues(eigenValues,indices);
    if(calcVectors){
      eigenVectors.Resize(n, nev);
      for (auto i = 0; i < nev; i++) {
        auto li = indices[i];
        for (auto x = 0; x < n; x++)
          eigenVectors(x, i) = eigenVectorsLapack(x, li);
      }
    }
  }else{
    if(calcVectors) eigenVectors = std::move(eigenVectorsLapack);
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

  TPZFMatrix <CTVar> eigenVectorsLapack;
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
  if(calcVectors) eigenVectorsLapack.Redim(n, n);
  
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
          eigenVectorsLapack( iVec , iCol) = z(iVec,iCol);
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
             (varfloatcomplex *)&eigenVectorsLapack(0,0), &ldz,
             (varfloatcomplex *)work.begin(),rwork.begin(), &info);
    }else if constexpr (std::is_same_v<TVar,std::complex<double>>){
      zhbgv_(&jobz, &uplo, &n, &ka, &kb,
             (vardoublecomplex *)A.fDiag.begin(), &ldab,
             (vardoublecomplex *)B.fDiag.begin(), &ldbb, w.begin(),
             (vardoublecomplex *)&eigenVectorsLapack(0,0), &ldz,
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
  const auto nev = this->NEigenpairs();
  if(nev > 0){
    TPZManVector<int,20> indices;
    this->SortEigenvalues(eigenValues,indices);
    eigenVectors.Resize(n,nev);
    for(auto i = 0; i < nev; i++){
      auto li = indices[i];
      for(auto x = 0; x < n; x++)
        eigenVectors(x,i) = eigenVectorsLapack(x,li); 
    }
  }else{
    eigenVectors = std::move(eigenVectorsLapack);
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