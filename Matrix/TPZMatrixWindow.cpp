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
void
TPZMatrixWindow<TVar>::MultAdd(const TPZMatrixWindow<TVar> &x,const TPZMatrixWindow<TVar> &y, TPZMatrixWindow<TVar> &z,
                               const TVar alpha,const TVar beta,const int opt_a, const int opt_x) const
{
  if(opt_x==0){
    //default checks, same as TPZMatrix<T>::MultAddChecks
    if ((!opt_a && this->Cols() != x.Rows()) || (opt_a && this->Rows() != x.Rows())) {
      TPZMatrix<TVar>::Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
      return;
    }
    if(beta != (TVar)0. && ((!opt_a && this->Rows() != y.Rows()) || (opt_a && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
      TPZMatrix<TVar>::Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
      return;
    }
    if(!opt_a) {
      if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
        z.Redim(this->Rows(),x.Cols());
      }
    } else {
      if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
        z.Redim(this->Cols(),x.Cols());
      }
    }
  }else{
    //checks for transposed x
    if ((!opt_a && this->Cols() != x.Cols()) || (opt_a && this->Rows() != x.Cols())) {
      TPZMatrix<TVar>::Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
      return;
    }
    if(beta != (TVar)0. && ((!opt_a && this->Rows() != y.Rows()) || (opt_a && this->Cols() != y.Rows()) || y.Cols() != x.Rows())) {
      TPZMatrix<TVar>::Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
      return;
    }
    if(!opt_a) {
      if(z.Cols() != x.Rows() || z.Rows() != this->Rows()) {
        z.Redim(this->Rows(),x.Rows());
      }
    } else {
      if(z.Cols() != x.Rows() || z.Rows() != this->Cols()) {
        z.Redim(this->Cols(),x.Rows());
      }
    }
  }
  DebugStop();
  // if (beta != (TVar)0) {
  //   z = y;
  // }else{
  //   z.Zero();
  // }

  const int64_t rows = this->Rows();
  const int64_t cols = this->Cols();
  const int64_t xrows = x.Rows();
  const int64_t xcols = x.Cols();
    
    
  if(rows == 0 || cols == 0 || xrows == 0 || xcols == 0){
    if (beta != (TVar)0) {
      z *= beta;
    }
    return;
  }

#ifdef USING_LAPACK
    
  //0: no transpose, 1: transpose, 2: conj trans
  const CBLAS_TRANSPOSE transp_a =
    opt_a == 0 ? CblasNoTrans : (opt_a == 1 ? CblasTrans : CblasConjTrans);
  const auto dim1 = opt_a == 0 ? rows : cols;
  const auto dim2 = opt_a == 0 ? cols : rows;

  const CBLAS_TRANSPOSE transp_x =
    opt_x == 0 ? CblasNoTrans : (opt_x == 1 ? CblasTrans : CblasConjTrans);
  
  if constexpr (std::is_same_v<TVar,double>){
    cblas_dgemm(CblasColMajor, transp_a, transp_x, dim1, xcols, dim2,
                alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, beta, z.fStorage, z.fLeadingDim);
    return;
  } else if constexpr (std::is_same_v<TVar,float>){
    cblas_sgemm(CblasColMajor, transp_a, transp_x, dim1, xcols, dim2,
                alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, beta, z.fStorage, z.fLeadingDim);
    return;
  } else if constexpr (std::is_same_v<TVar,std::complex<double>>){
    cblas_zgemm(CblasColMajor, transp_a, transp_x, dim1, xcols, dim2,
                &alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, &beta, z.fStorage, z.fLeadingDim);
    return;
  } else if constexpr (std::is_same_v<TVar,std::complex<float>>){
    cblas_cgemm(CblasColMajor, transp_a, transp_x, dim1, xcols, dim2,
                &alpha, this->fStorage, this->fLeadingDim, x.fStorage,
                x.fLeadingDim, &beta, z.fStorage, z.fLeadingDim);
    return;
  }
#endif
  DebugStop();
  // if (beta != (TVar)0) {
  //   z *= beta;
  // }
  // for (auto ic = 0; ic < xcols; ic++) {
  //   if(!opt) {
  //     for (auto c = 0; c<cols; c++) {
  //       TVar * zp = &z(0,ic), *zlast = zp+rows;
  //       TVar * fp = fElem +rows*c;
  //       const TVar * xp = &x.g(c,ic);
  //       while(zp < zlast) {
  //         *zp += alpha* *fp++ * *xp;
  //         zp ++;
  //       }
  //     }
  //   } else {
  //     TVar * fp = fElem,  *zp = &z(0,ic);
  //     for (auto c = 0; c<cols; c++) {
  //       TVar val = 0.;
  //       // bug correction philippe 5/2/97
  //       //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
  //       const TVar *xp = &x.g(0,ic);
  //       const TVar *xlast = xp + rows;
  //       if constexpr (is_complex<TVar>::value){
  //         if(opt==2){
  //           while(xp < xlast) {
  //             val += std::conj(*fp++) * *xp++;
  //           }
  //           *zp++ += alpha *val;
  //           continue;//continue from the for loop
  //         }
  //       }
  //       while(xp < xlast) {
  //         val += *fp++ * *xp++;
  //       }
  //       *zp++ += alpha *val;
  //     }
  //   }
  // }
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