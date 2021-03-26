/**
 * @file
 */

#include "TExtFunction.h"

template<class TVar>
TExtFunction<TVar>::TExtFunction() {

}

template<class TVar>
TExtFunction<TVar>::~TExtFunction() {

}
    
template<class TVar>
void TExtFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
  f.Resize(1);
  df.Resize(2,1);
  double r = 2.-x[0];
  if(r < 1e-6) r = 1e-6;
  f[0] = (TVar)log(r);
  df(0,0) = (TVar)(-1./(r));
  df(1,0) = (TVar)0.;
}  

template<class TVar>
int TExtFunction<TVar>::NFunctions() const {
  return 1;
}
  
template<class TVar>
int TExtFunction<TVar>::PolynomialOrder() const{
  return 1;
}

template class TExtFunction<float>;
template class TExtFunction<double>;
template class TExtFunction<long double>;