//$Id: TDiscoFunction.cpp,v 1.3 2008-04-04 13:36:47 fortiago Exp $

#include "TDiscoFunction.h"

 TDiscoFunction:: TDiscoFunction(){

}

 TDiscoFunction::~ TDiscoFunction(){

}
    
void  TDiscoFunction::Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df){
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

int  TDiscoFunction::NFunctions(){
  return 1;
}
  
int  TDiscoFunction::PolynomialOrder(){
  return 100;
}

