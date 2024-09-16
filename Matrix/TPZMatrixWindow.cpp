#include "TPZMatrixWindow.h"
#include "pzfmatrix.h"

#ifdef USING_LAPACK
#include "TPZLapack.h"
#endif

template<class TVar>
TPZMatrixWindow<TVar>::TPZMatrixWindow(TVar* mem_area, const int nrows, const int ncols,
                                       const int leading_dim, const int size_mem)
  : fStorage(mem_area), fLeadingDim(leading_dim)
{
  const int nr_orig = leading_dim;
  const int nc_orig = size_mem/leading_dim;
  CheckConstructor(0,0,nrows,ncols,nr_orig,nc_orig);
  this->fRow = nrows;
  this->fCol = ncols;
}

template<class TVar>
TPZMatrixWindow<TVar>::TPZMatrixWindow(TPZFMatrix<TVar> &mat, const int i, const int j,
                                       const int nrows, const int ncols)
  : fLeadingDim(mat.Rows())
{
  CheckConstructor(i,j,nrows,ncols,mat.Rows(),mat.Cols());
  //pointer of first window position
  fStorage = mat.Elem() + mat.Rows()*j+i;
  this->fRow = nrows;
  this->fCol = ncols;
}

template<class TVar>
TPZMatrixWindow<TVar> &TPZMatrixWindow<TVar>::operator=(const TPZFMatrix<TVar> &mat){
  if(this->Rows()!=mat.Rows() || this->Cols() != mat.Cols()){
    DebugStop();
  }
  const auto m_leading_dim = mat.Rows();
  const auto nc = this->Cols();
  const auto nr = this->Rows();
  //mem dist to advance between cols
  const auto my_adv = this->fLeadingDim;
  const auto their_adv = m_leading_dim;
  
  int64_t my_count{0},their_count{0};
  for(int ic = 0; ic < nc; ic++){
    auto *my_begin = &this->fStorage[my_count];
    auto *their_begin = &mat.Elem()[their_count];
    auto *their_end = their_begin+nr;
    std::copy(their_begin,their_end,my_begin);
    my_count+=my_adv;
    their_count+=their_adv;
  }
  // const auto sz = nc*nr;
  // for(int i = 0; i < sz; i++){
  //   const int ic = i/nr;
  //   const int ir = i%nr;
  //   this->fStorage[ic*this->fLeadingDim+ir] = mat.Elem()[ic*m_leading_dim+ir];
  // }
  
  return *this;
}
template<class TVar>
TPZMatrixWindow<TVar> &TPZMatrixWindow<TVar>::operator=(const TPZMatrixWindow<TVar> &mat){
  if(this->Rows()!=mat.Rows() || this->Cols() != mat.Cols()){
    DebugStop();
  }
  const auto m_leading_dim = mat.fLeadingDim;
  const auto nc = this->Cols();
  const auto nr = this->Rows();
  //mem dist to advance between cols
  const auto my_adv = this->fLeadingDim;
  const auto their_adv = m_leading_dim;

  int64_t my_count{0},their_count{0};
  for(int ic = 0; ic < nc; ic++){
    auto *my_begin = &this->fStorage[my_count];
    auto *their_begin = &mat.fStorage[their_count];
    auto *their_end = their_begin+nr;
    std::copy(their_begin,their_end,my_begin);
    my_count+=my_adv;
    their_count+=their_adv;
  }
  // const auto sz = nc*nr;
  // for(int i = 0; i < sz; i++){
  //   const int ic = i/nr;
  //   const int ir = i%nr;
  //   this->fStorage[ic*this->fLeadingDim+ir] = mat.fStorage[ic*m_leading_dim+ir];
  // }
  return *this;
}

template<class TVar>
TPZMatrixWindow<TVar> &TPZMatrixWindow<TVar>::operator*=(const TVar val){
  const auto nc = this->Cols();
  const auto nr = this->Rows();
  auto *my_ptr = &this->fStorage[0];
  //mem dist to advance between cols
  const auto my_adv = this->fLeadingDim-nr;
  
  for(int ic = 0; ic < nc; ic++){
    for(auto ir = 0; ir < nr; ir++,my_ptr++){
      *my_ptr *= val;
    }
    my_ptr+=my_adv;
  }
  // const auto sz = nc*nr;
  // for(int i = 0; i < sz; i++){
  //   const int ic = i/nr;
  //   const int ir = i%nr;
  //   this->fStorage[ic*this->fLeadingDim+ir] *= val;
  // }
  return *this;
}

// implement the matrix vector product z = alpha * opt(this)*opt(x) + beta * y
template<class TVar>
void
TPZMatrixWindow<TVar>::MultAdd(const TPZMatrixWindow<TVar> &x,const TPZMatrixWindow<TVar> &y, TPZMatrixWindow<TVar> &z,
                               const TVar alpha,const TVar beta,const int opt_a, const int opt_x) const
{
  if(opt_x==0){
    //default checks, same as TPZMatrix<T>::MultAddChecks
    if ((!opt_a && this->Cols() != x.Rows()) || (opt_a && this->Rows() != x.Rows())) {
      PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix x with incompatible dimensions"<<std::endl;
      DebugStop();
    }
    if(beta != (TVar)0. && ((!opt_a && this->Rows() != y.Rows()) || (opt_a && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
      PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix y with incompatible dimensions"<<std::endl;
      DebugStop();
    }
    if(!opt_a) {
      if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
        PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix z with incompatible dimensions"<<std::endl;
      DebugStop();
      }
    } else {
      if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
        PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix z with incompatible dimensions"<<std::endl;
      DebugStop();
      }
    }
  }else{
    //checks for transposed x
    if ((!opt_a && this->Cols() != x.Cols()) || (opt_a && this->Rows() != x.Cols())) {
      PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix x with incompatible dimensions"<<std::endl;
      DebugStop();      
    }
    if(beta != (TVar)0. && ((!opt_a && this->Rows() != y.Rows()) || (opt_a && this->Cols() != y.Rows()) || y.Cols() != x.Rows())) {
      PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix y with incompatible dimensions"<<std::endl;
      DebugStop();
    }
    if(!opt_a) {
      if(z.Cols() != x.Rows() || z.Rows() != this->Rows()) {
        PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix z with incompatible dimensions"<<std::endl;
      DebugStop();
      }
    } else {
      if(z.Cols() != x.Rows() || z.Rows() != this->Cols()) {
        PZError<<__PRETTY_FUNCTION__
             <<"\n:Matrix z with incompatible dimensions"<<std::endl;
      DebugStop();
      }
    }
  }

  if (beta != (TVar)0) {
    z = y;
  }else{
    z.Zero();
  }

  const int64_t rows = this->Rows();
  const int64_t cols = this->Cols();
  const int64_t xrows = x.Rows();
  const int64_t xcols = x.Cols();

  ////we do not allow zero-sized windows - why?
  // if(rows == 0 || cols == 0 || xrows == 0 || xcols == 0){
  //   if (beta != (TVar)0) {
  //     z *= beta;
  //   }
  //   return;
  // }

#ifdef USING_LAPACK
    
  //0: no transpose, 1: transpose, 2: conj trans
  const CBLAS_TRANSPOSE transp_a =
    opt_a == 0 ? CblasNoTrans : (opt_a == 1 ? CblasTrans : CblasConjTrans);
  const auto dim1 = opt_a == 0 ? rows : cols;
  const auto dim2 = opt_a == 0 ? cols : rows;
  const auto dim3 = opt_x == 0 ? xcols : xrows;
  const CBLAS_TRANSPOSE transp_x =
    opt_x == 0 ? CblasNoTrans : (opt_x == 1 ? CblasTrans : CblasConjTrans);
  
  if constexpr (std::is_same_v<TVar,double>){
    cblas_dgemm(CblasColMajor, transp_a, transp_x, dim1, dim3, dim2,
                alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, beta, z.fStorage, z.fLeadingDim);
    return;
  } else if constexpr (std::is_same_v<TVar,float>){
    cblas_sgemm(CblasColMajor, transp_a, transp_x, dim1, dim3, dim2,
                alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, beta, z.fStorage, z.fLeadingDim);
    return;
  } else if constexpr (std::is_same_v<TVar,std::complex<double>>){
    cblas_zgemm(CblasColMajor, transp_a, transp_x, dim1, dim3, dim2,
                &alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, &beta, z.fStorage, z.fLeadingDim);
    return;
  } else if constexpr (std::is_same_v<TVar,std::complex<float>>){
    cblas_cgemm(CblasColMajor, transp_a, transp_x, dim1, dim3, dim2,
                &alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, &beta, z.fStorage, z.fLeadingDim);
    return;
  }
#endif
  if (beta != (TVar)0) {
    z *= beta;
  }
  if(!opt_x) {
    for (auto ic = 0; ic < xcols; ic++) {
      if(!opt_a) {
        for (auto c = 0; c<cols; c++) {
          TVar * zp = &z(0,ic), *zlast = zp+rows;
          TVar * fp = &g(0,c);//fElem +rows*c;
          const TVar * xp = &x.g(c,ic);
          while(zp < zlast) {
            *zp += alpha* *fp++ * *xp;
            zp ++;
          }
        }
      } else if(opt_a) {
        for (auto c = 0; c<cols; c++) {
          TVar *zp = &z(c,ic);
          TVar * fp = &g(0,c);
          TVar val = 0.;
          // bug correction philippe 5/2/97
          //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
          const TVar *xp = &x.g(0,ic);
          const TVar *xlast = xp + xrows;
          if constexpr (is_complex<TVar>::value){
            if(opt_a==2){
              while(xp < xlast) {
                val += std::conj(*fp++) * *xp++;
              }
              *zp += alpha *val;
              continue;//continue from the for loop
            }
            else if(opt_a==1){
              while(xp < xlast) {
                val += *fp++ * *xp++;
              }
              *zp += alpha *val;
              continue;//continue from the for loop
            }
          } else {
            while(xp < xlast) {
              val += *fp++ * *xp++;
            }
            *zp += alpha *val;
          }
        }
      }
    }
  } else if (opt_x) {
    if(opt_a) {
      for (auto xr = 0; xr < xrows; xr++) {
        for (auto c = 0; c<cols; c++) {
          TVar * zp = &z(c,xr);
          TVar * fp = &g(0,c);//fElem +rows*c;
          TVar * fplast = fp + rows;
          const TVar * xp = &x.g(xr,0);
          if constexpr (is_complex<TVar>::value){
            if(opt_a==2 && opt_x == 1){
              TVar val = 0.;
              while(fp < fplast) {
                val += alpha* std::conj(*fp++) * *xp;
                xp += x.fLeadingDim;
              }
              *zp += val;
            } else if(opt_a==2 && opt_x == 2){
              TVar val = 0.;
              while(fp < fplast) {
                val += alpha* std::conj(*fp++) * std::conj(*xp);
                xp += x.fLeadingDim;
              }
              *zp += val;
            } else if(opt_a == 1 && opt_x ==2) {
              TVar val = 0.;
              while(fp < fplast) {
                val += alpha* *fp++ * std::conj(*xp);
                xp += x.fLeadingDim;
              }
              *zp += val;
            } else {
              TVar val = 0.;
              while(fp < fplast) {
                val += alpha* *fp++ * *xp;
                xp += x.fLeadingDim;
              }
              *zp += val;
            }
          } else {
            TVar val = 0.;
            while(fp < fplast) {
              val += alpha* *fp++ * *xp;
              xp += x.fLeadingDim;
            }
            *zp += val;
          }
        }
      } 
    } else if(opt_a == 0) {
      for (auto xr = 0; xr<x.Rows(); xr++) {
        for (auto c = 0; c < cols; c++) {
          TVar * fp = &g(0,c);
          TVar * fplast = fp + rows;
          TVar *zp = &z(0,xr);
          // bug correction philippe 5/2/97
          //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
          // const TVar *xlast = xp + x.fLeadingDim;
          if constexpr (is_complex<TVar>::value){
            if(opt_x==2){
              TVar val = alpha *std::conj(x.g(xr,c));
              while(fp < fplast) {
                *zp++ += *fp++ * val;
              }
              continue;//continue from the for loop
            } else if (opt_x==1){
              TVar val = alpha *x.g(xr,c);
              while(fp < fplast) {
                *zp++ += *fp++ * val;
              }
              continue;//continue from the for loop
            }
          } else {
            TVar val = alpha*x.g(xr,c);
            while(fp < fplast) {
              *zp++ += *fp++ * val;
            }
          }
        }
      }
    }
  }
}

template<class TVar>
void
TPZMatrixWindow<TVar>::CheckConstructor(const int i, const int j, const int nr, const int nc,
                                        const int nr_orig, const int nc_orig){
  if(nr < 1 || nc < 1 || i+nr>nr_orig || j+nc>nc_orig || i < 0 || j < 0){
    PZError<<__PRETTY_FUNCTION__
           <<"\nInvalid parameters:\n"
           <<"i "<<i<<'\n'
           <<"j "<<j<<'\n'
           <<"nr "<<nr<<'\n'
           <<"nc "<<nc<<'\n'
           <<"nr_orig "<<nr_orig<<'\n'
           <<"nc_orig "<<nc_orig<<std::endl;
    DebugStop();
  }
}

template class TPZMatrixWindow<float >;
template class TPZMatrixWindow<double >;
template class TPZMatrixWindow<long double>;

template class TPZMatrixWindow< std::complex<float> >;
template class TPZMatrixWindow< std::complex<double> >;
template class TPZMatrixWindow< std::complex<long double> >;