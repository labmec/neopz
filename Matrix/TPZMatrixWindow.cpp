#include "TPZMatrixWindow.h"
#include "pzfmatrix.h"

template<class TVar>
TPZMatrixWindow<TVar>::TPZMatrixWindow(TVar* mem_area, const int nrows, const int ncols, const int leading_dim)
  : fStorage(mem_area), fLeadingDim(leading_dim)
{
  this->fRow = nrows;
  this->fCol = ncols;
}

template<class TVar>
TPZMatrixWindow<TVar>::TPZMatrixWindow(TPZFMatrix<TVar> &mat, const int i, const int j,
                                       const int nrows, const int ncols)
  : fLeadingDim(mat.Rows())
{
  //pointer of first window position
  fStorage = mat.Elem() + mat.Rows()*j+i;
  this->fRow = nrows;
  this->fCol = ncols;
}

template<class TVar>
void
TPZMatrixWindow<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                               const TVar alpha,const TVar beta,const int opt_a, const int opt_x) const
{
  //implement me
  DebugStop();
}


template class TPZMatrixWindow<float >;
template class TPZMatrixWindow<double >;
template class TPZMatrixWindow<long double>;

template class TPZMatrixWindow< std::complex<float> >;
template class TPZMatrixWindow< std::complex<double> >;
template class TPZMatrixWindow< std::complex<long double> >;