//$Id: TDiscoFunction.cpp,v 1.3 2008-04-04 13:36:47 fortiago Exp $

#include "TDiscoFunction.h"

template<class TVar>
TDiscoFunction<TVar>::TDiscoFunction(){

}

template<class TVar>
TDiscoFunction<TVar>::~TDiscoFunction(){

}
    
template<class TVar>
void  TDiscoFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
  f.Resize(1);
  df.Resize(8,1);
  df.Zero();
  double Xc = 0.;
  double Yc = 0.;
  double Xp = x[0];
  double Yp = x[1];
  double r = sqrt( (Xp-Xc)*(Xp-Xc) +(Yp-Yc)*(Yp-Yc) );
  f[0] = log(r);
  df(0,0) = (Xp - Xc)/(r*r);
  df(1,0) = (Yp - Yc)/(r*r);
  df(2,0) = 0.; //Laplaciano
  ///as outras derivadas nao serao usadas
}  

template<class TVar>
int  TDiscoFunction<TVar>::NFunctions()const {
  return 1;
}
  
template<class TVar>
int  TDiscoFunction<TVar>::PolynomialOrder() const{
  return 100;
}

template class TDiscoFunction<float>;
template class TDiscoFunction<double>;
template class TDiscoFunction<long double>;